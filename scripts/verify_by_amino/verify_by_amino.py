import argparse
import os
import sys
import time
from pathlib import Path

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

STOP_SYMBOL = "*"
MATCH_IDENTITY_THRESHOLD = 0.5

DEFAULT_AA_REFERENCE = str(
    Path(__file__).resolve().parent.parent.parent / "data"
    / "HomoSapiens_IMGTGENEDB-ReferenceSequences.fasta"
)


def read_fasta(path):
    """Читает fasta-файл, сохраняя заголовки как есть (без разбора по '|')."""
    records = []
    header = None
    chunks = []

    def flush():
        if header is not None:
            records.append((header, "".join(chunks)))

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line[1:]
                chunks = []
            else:
                chunks.append(line.strip())
        flush()
    return records


def write_fasta(records, path):
    with open(path, "w", encoding="utf-8") as f:
        for header, seq in records:
            f.write(f">{header}\n{seq}\n")


def trim_to_complete_codons(nt_seq, frame_offset):
    seq = nt_seq[frame_offset:]
    complete_len = len(seq) - (len(seq) % 3)
    return seq[:complete_len]


def translate_frame(nt_seq, frame_offset):
    """Транслирует нуклеотидную последовательность в заданной рамке считывания.

    Возвращает (protein, trimmed_nt) или (None, None), если после обрезки
    не осталось ни одного целого кодона. trimmed_nt — это nt_seq, обрезанная
    под рамку: без сдвига в начале и без "хвостовых" 1-2 нуклеотидов в конце
    (то есть длина trimmed_nt всегда кратна 3).
    """
    trimmed = trim_to_complete_codons(nt_seq, frame_offset)
    trimmed = trimmed.replace("-", "")
    if len(trimmed) < 3:
        return None, None
    protein = str(Seq(trimmed).translate(to_stop=False))
    return protein, trimmed


def classify_stop_codon(protein):
    stop_positions = [i for i, aa in enumerate(protein) if aa == STOP_SYMBOL]
    if not stop_positions:
        return "no_stop"
    if len(stop_positions) == 1 and stop_positions[0] == len(protein) - 1:
        return "stop_at_end"
    return "premature_stop"


def strip_trailing_stop(protein):
    return protein[:-1] if protein.endswith(STOP_SYMBOL) else protein


def read_aa_reference(path):
    """Читает AA-справочник IMGT."""
    seqs = {}
    header = None
    chunks = []
    total_records = 0

    def flush():
        nonlocal total_records
        if header is not None:
            total_records += 1
            parts = header.split("|")
            name = parts[1].strip() if len(parts) > 1 else header.strip()
            seqs[name] = "".join(chunks).upper()

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
        flush()
    return seqs, total_records


def best_reference_match(protein, reference_dict):
    """Ищет наиболее похожую последовательность в справочнике."""
    aligner = PairwiseAligner()
    aligner.mode = "local"

    best_name = None
    best_score = float("-inf")
    best_ref_len = 0
    for name, ref_seq in reference_dict.items():
        score = aligner.score(protein, ref_seq)
        if score > best_score:
            best_score = score
            best_name = name
            best_ref_len = len(ref_seq)
    return best_name, best_score, best_ref_len


def is_confident_match(best_score, best_ref_len, threshold=MATCH_IDENTITY_THRESHOLD):
    if not best_ref_len:
        return False
    return (best_score / best_ref_len) >= threshold


def process_sequence(nt_seq, reference_dict, counters):
    """Пробует рамки считывания по очереди (1, 2, 3).

    Возвращает (status, protein, trimmed_nt, frame_number):
        ("ok", protein, обрезанная_под_рамку_ДНК, номер_рамки)
        ("no_match", None, None, None)
        ("premature_stop", None, None, None)
    """
    had_valid_stop_frame = False

    for frame_offset in (0, 1, 2):
        protein, trimmed_nt = translate_frame(nt_seq, frame_offset)
        if protein is None:
            continue

        stop_status = classify_stop_codon(protein)
        if stop_status == "premature_stop":
            continue

        had_valid_stop_frame = True
        clean_protein = strip_trailing_stop(protein)

        if not reference_dict:
            counters[stop_status] += 1
            counters["matched"] += 1
            return "ok", clean_protein, trimmed_nt, frame_offset + 1

        _, best_score, best_ref_len = best_reference_match(clean_protein, reference_dict)
        if is_confident_match(best_score, best_ref_len):
            counters[stop_status] += 1
            counters["matched"] += 1
            return "ok", clean_protein, trimmed_nt, frame_offset + 1

    if had_valid_stop_frame:
        counters["no_match"] += 1
        return "no_match", None, None, None

    counters["premature_stop"] += 1
    return "premature_stop", None, None, None


def process_file(input_path, output_dir, quarantine_dir, no_match_dir, reference_dict, counters):
    """Обрабатывает один fasta-файл. На выход — нуклеотиды прошедших проверку."""
    records = read_fasta(input_path)
    filename = os.path.basename(input_path)
    print(f"\n--- {filename}: {len(records)} последовательностей ---")

    ok_records = []
    premature_records = []
    no_match_records = []

    for i, (header, nt_seq) in enumerate(records, 1):
        nt_seq = nt_seq.upper()
        status, protein, trimmed_nt, _frame = process_sequence(nt_seq, reference_dict, counters)

        if status == "ok":
            ok_records.append((header, trimmed_nt))
        elif status == "no_match":
            no_match_records.append((header, nt_seq))
        else:
            premature_records.append((header, nt_seq))

        if i % 20 == 0 or i == len(records):
            print(f"  обработано {i}/{len(records)}")

    if ok_records:
        write_fasta(ok_records, os.path.join(output_dir, filename))
    if premature_records:
        write_fasta(premature_records, os.path.join(quarantine_dir, filename))
    if no_match_records:
        write_fasta(no_match_records, os.path.join(no_match_dir, filename))

    print(f"  ok={len(ok_records)}, no_match={len(no_match_records)}, premature_stop={len(premature_records)}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Верификация нуклеотидных fasta через трансляцию и сверку по IMGT. На выход — нуклеотиды."
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Входная папка с nt-fasta файлами (vj_filtered/)")
    parser.add_argument("-o", "--output", required=True,
                        help="Выходная папка для прошедших проверку nt-fasta")
    parser.add_argument("--aa-reference", default=DEFAULT_AA_REFERENCE,
                        help="Путь к AA-справочнику IMGT без гэпов")
    return parser.parse_args()


def main():
    args = parse_args()

    input_dir = args.input
    output_dir = args.output
    quarantine_dir = os.path.join(output_dir, "premature_stop")
    no_match_dir = os.path.join(output_dir, "no_match")

    if not os.path.isdir(input_dir):
        print(f"Ошибка: папка {input_dir} не найдена.", file=sys.stderr)
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(quarantine_dir, exist_ok=True)
    os.makedirs(no_match_dir, exist_ok=True)

    reference_path = args.aa_reference
    if os.path.isfile(reference_path):
        reference_dict, total_ref = read_aa_reference(reference_path)
        print(f"Загружен AA-справочник: {len(reference_dict)} уникальных генов, {total_ref} записей")
    else:
        reference_dict = {}
        print(
            f"Внимание: AA-справочник не найден по пути {reference_path}. "
            "Сверка со справочником пропускается — берём первый вариант, "
            "прошедший проверку по стоп-кодону."
        )

    counters = {
        "no_stop": 0,
        "stop_at_end": 0,
        "matched": 0,
        "no_match": 0,
        "premature_stop": 0,
    }

    fasta_files = sorted(f for f in os.listdir(input_dir) if f.lower().endswith((".fasta", ".fa")))
    print(f"Найдено fasta-файлов: {len(fasta_files)}")
    t_start = time.time()

    for idx, filename in enumerate(fasta_files, 1):
        print(f"\n[{idx}/{len(fasta_files)}] {filename}")
        process_file(
            os.path.join(input_dir, filename),
            output_dir, quarantine_dir, no_match_dir,
            reference_dict, counters,
        )
        elapsed = time.time() - t_start
        print(f"  прошло {elapsed:.1f} сек")

    elapsed = time.time() - t_start
    print(f"\nГотово! ({elapsed:.1f} сек)")
    print(f"Без стоп-кодона: {counters['no_stop']}")
    print(f"Стоп-кодон в конце: {counters['stop_at_end']}")
    print(f"Совпало со справочником (сохранено в verify_by_amino/): {counters['matched']}")
    print(f"Не нашли совпадения ни в одной рамке (в no_match/): {counters['no_match']}")
    print(f"Стоп-кодон не в конце во всех рамках (в premature_stop/): {counters['premature_stop']}")


if __name__ == "__main__":
    main()