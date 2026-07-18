import argparse
import os
import sys
import time

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

# ============ ПУТИ ============
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# paths.py лежит на уровень выше (в scripts/), а не рядом с этим файлом
SCRIPTS_DIR = os.path.dirname(BASE_DIR)
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from paths import get_paths  # noqa: E402

STOP_SYMBOL = "*"
MATCH_IDENTITY_THRESHOLD = 0.5  # доля совпавших позиций от длины белка, чтобы считать совпадением

DEFAULT_INPUT_SUBFOLDER = "vj_filtered"
DEFAULT_OUTPUT_SUBFOLDER = "amino"
DEFAULT_REFERENCE_FILENAME = "HomoSapiens_IMGTGENEDB-ReferenceSequences.fasta"


def read_fasta(path):
    """Читает fasta-файл, сохраняя заголовки как есть (без разбора по '|') —
    нам нужно сохранить исходный заголовок 1-в-1 при трансляции."""
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
    """Обрезает последовательность под рамку считывания: убирает frame_offset
    нуклеотидов слева и 1-2 "хвостовых" нуклеотида справа, если они не
    складываются в целый кодон.
    """
    seq = nt_seq[frame_offset:]
    complete_len = len(seq) - (len(seq) % 3)
    return seq[:complete_len]


def translate_frame(nt_seq, frame_offset):
    """Транслирует нуклеотидную последовательность в заданной рамке считывания
    (0, 1 или 2). Возвращает None, если после обрезки не осталось ни одного
    целого кодона.
    """
    trimmed = trim_to_complete_codons(nt_seq, frame_offset)
    trimmed = trimmed.replace("-", "")
    if len(trimmed) < 3:
        return None
    return str(Seq(trimmed).translate(to_stop=False))


def classify_stop_codon(protein):
    """Определяет, где встречается стоп-кодон в транслированном белке.

    Возвращает:
        'no_stop'         — стоп-кодона нет вообще
        'stop_at_end'      — стоп-кодон есть, ровно один, и он в самом конце
        'premature_stop'   — стоп-кодон встретился раньше конца (или их несколько)
    """
    stop_positions = [i for i, aa in enumerate(protein) if aa == STOP_SYMBOL]
    if not stop_positions:
        return "no_stop"
    if len(stop_positions) == 1 and stop_positions[0] == len(protein) - 1:
        return "stop_at_end"
    return "premature_stop"


def strip_trailing_stop(protein):
    """Убирает завершающий '*', если он есть."""
    return protein[:-1] if protein.endswith(STOP_SYMBOL) else protein


def read_aa_reference(path):
    """Читает AA-справочник IMGT: ключ — имя гена/аллели из 2-го поля
    заголовка (между '|'), значение — аминокислотная последовательность."""
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
    """Ищет наиболее похожую последовательность в справочнике методом
    локального попарного выравнивания.

    Возвращает (имя_гена, скор, длина_гена_из_справочника). Длина гена нужна
    отдельно, потому что наш транслированный белок (V+junction+J) длиннее
    любого одиночного V- или J-гена из справочника — сравнивать долю
    совпадения нужно относительно длины СПРАВОЧНОГО гена, а не всего белка.
    """
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
    """Скор локального выравнивания (при match=1) примерно равен числу
    совпавших позиций. Делим на длину гена из справочника (не на длину всего
    транслированного белка), чтобы получить долю покрытия этого гена."""
    if not best_ref_len:
        return False
    return (best_score / best_ref_len) >= threshold


def process_sequence(nt_seq, reference_dict, counters):
    """Пробует рамки считывания по очереди (1, 2, 3).

    Возвращает (status, protein, frame_number):
        ("ok", protein, номер_рамки)   — нашли подходящий вариант
        ("no_match", None, None)        — стоп-кодон был валиден хоть в одной
                                           рамке, но ни одна не совпала со
                                           справочником
        ("premature_stop", None, None)  — ни в одной рамке стоп-кодон не
                                           оказался в конце (везде преждевременный)
    """
    had_valid_stop_frame = False

    for frame_offset in (0, 1, 2):
        protein = translate_frame(nt_seq, frame_offset)
        if protein is None:
            continue  # нечего транслировать в этой рамке

        stop_status = classify_stop_codon(protein)
        if stop_status == "premature_stop":
            continue  # эта рамка не годится, пробуем следующую

        had_valid_stop_frame = True
        clean_protein = strip_trailing_stop(protein)

        if not reference_dict:
            # справочника нет — берём первый вариант, прошедший проверку по стоп-кодону
            counters[stop_status] += 1
            counters["matched"] += 1
            return "ok", clean_protein, frame_offset + 1

        _, best_score, best_ref_len = best_reference_match(clean_protein, reference_dict)
        if is_confident_match(best_score, best_ref_len):
            counters[stop_status] += 1
            counters["matched"] += 1
            return "ok", clean_protein, frame_offset + 1

    if had_valid_stop_frame:
        counters["no_match"] += 1
        return "no_match", None, None

    counters["premature_stop"] += 1
    return "premature_stop", None, None


def process_file(input_path, output_dir, quarantine_dir, no_match_dir, reference_dict, counters):
    """Обрабатывает один fasta-файл, раскладывая записи по трём возможным
    итогам: успешно транслированные, преждевременный стоп-кодон, нет совпадения."""
    records = read_fasta(input_path)
    filename = os.path.basename(input_path)
    print(f"\n--- {filename}: {len(records)} последовательностей ---")

    ok_records = []
    premature_records = []
    no_match_records = []

    for i, (header, nt_seq) in enumerate(records, 1):
        nt_seq = nt_seq.upper()
        status, protein, _frame = process_sequence(nt_seq, reference_dict, counters)

        if status == "ok":
            ok_records.append((header, protein))
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
        description="Трансляция нуклеотидных fasta (V+J) в аминокислотные со сверкой по IMGT."
    )
    parser.add_argument(
        "-k", "--key",
        default="BCR",
        help="Название подпапки внутри results/ (по умолчанию: BCR)",
    )
    parser.add_argument(
        "--input-subfolder",
        default=DEFAULT_INPUT_SUBFOLDER,
        help=f"Подпапка внутри results/<key>/ с входными nt-fasta (по умолчанию: {DEFAULT_INPUT_SUBFOLDER})",
    )
    parser.add_argument(
        "--output-subfolder",
        default=DEFAULT_OUTPUT_SUBFOLDER,
        help=f"Подпапка внутри results/<key>/ для результатов (по умолчанию: {DEFAULT_OUTPUT_SUBFOLDER})",
    )
    parser.add_argument(
        "--reference-fasta",
        default=None,
        help=(
            "Путь к AA-справочнику IMGT без гэпов. Если не указан, берётся "
            f"data/{DEFAULT_REFERENCE_FILENAME}"
        ),
    )
    return parser.parse_args()


def main():
    args = parse_args()
    paths = get_paths(args.key)

    input_dir = os.path.join(paths["output_dir"], args.input_subfolder)
    output_dir = os.path.join(paths["output_dir"], args.output_subfolder)
    quarantine_dir = os.path.join(output_dir, "premature_stop")
    no_match_dir = os.path.join(output_dir, "no_match")

    if not os.path.isdir(input_dir):
        print(f"Ошибка: путь {input_dir} не найден. Проверь --key и --input-subfolder.")
        return

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(quarantine_dir, exist_ok=True)
    os.makedirs(no_match_dir, exist_ok=True)


# !!!!!!!!!!!!!!!!!!!!
    reference_path = args.reference_fasta or os.path.join(os.path.dirname(paths["input_dir"]), DEFAULT_REFERENCE_FILENAME)
    print(args.reference_fasta)
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

    fasta_files = [f for f in os.listdir(input_dir) if f.lower().endswith((".fasta", ".fa"))]
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

    print("\nГотово!")
    print(f"Без стоп-кодона: {counters['no_stop']}")
    print(f"Стоп-кодон в конце: {counters['stop_at_end']}")
    print(f"Совпало со справочником (сохранено в {args.output_subfolder}/): {counters['matched']}")
    print(f"Не нашли совпадения ни в одной рамке (в no_match/): {counters['no_match']}")
    print(f"Стоп-кодон не в конце во всех рамках (в premature_stop/): {counters['premature_stop']}")


if __name__ == "__main__":
    main()
