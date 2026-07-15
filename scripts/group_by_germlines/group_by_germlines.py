import os
import re

import pandas as pd
from Bio import Align

# ============ ПУТИ ============
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# V/J справочники — свои для каждого локуса, D — только у тяжёлой цепи
# (у каппа/лямбда лёгких цепей нет D-сегмента, V соединяется с J напрямую).
LOCUS_VJ_FASTA = {
    "IGH": {"v": "IGHV.fasta", "j": "IGHJ.fasta"},
    "IGK": {"v": "IGKV.fasta", "j": "IGKJ.fasta"},
    "IGL": {"v": "IGLV.fasta", "j": "IGLJ.fasta"},
}
D_FASTA = "IGHD.fasta"

INPUT_FILE = os.path.join(BASE_DIR, "BCR_data_filtered.tsv")  # результат первого скрипта

OUT_ROOT = os.path.join(BASE_DIR, "grouped_by_germlines")

aligner = Align.PairwiseAligner()
aligner.mode = "local"  # локальное — ищем лучший совпадающий участок


def read_germline_fasta(path):
    """Читает IMGT fasta -> словарь {имя_гена: последовательность (в верхнем регистре)}."""
    seqs = {}
    name = None
    chunks = []

    def flush():
        if name is not None:
            seqs[name] = "".join(chunks).upper()

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                parts = line[1:].split("|")
                name = parts[1] if len(parts) > 1 else parts[0]
                chunks = []
            else:
                chunks.append(line)
        flush()
    return seqs


def safe_filename(name):
    return re.sub(r'[^A-Za-z0-9_\-\.]', '_', name)


def best_germline_match(query_seq, germline_dict):
    """Находит ген с максимальным скором выравнивания (наиболее похожий гермлайн)."""
    best_name = None
    best_score = float("-inf")
    query_seq = query_seq.upper()
    for name, gseq in germline_dict.items():
        score = aligner.score(query_seq, gseq)
        if score > best_score:
            best_score = score
            best_name = name
    return best_name, best_score


def build_folder(out_root, subfolder):
    path = os.path.join(out_root, subfolder)
    os.makedirs(path, exist_ok=True)
    return path


def write_fasta(records, out_dir, gene_name):
    fname = safe_filename(gene_name) + ".fasta"
    fpath = os.path.join(out_dir, fname)
    with open(fpath, "w", encoding="utf-8") as f:
        for seq_id, seq in records:
            f.write(f">{seq_id}\n{seq}\n")


def main():
    print("Читаю germline-справочники (полные последовательности)...")
    v_seqs_by_locus = {}
    j_seqs_by_locus = {}
    for locus, files in LOCUS_VJ_FASTA.items():
        v_seqs_by_locus[locus] = read_germline_fasta(os.path.join(BASE_DIR, files["v"]))
        j_seqs_by_locus[locus] = read_germline_fasta(os.path.join(BASE_DIR, files["j"]))
        print(f"  {locus}: V={len(v_seqs_by_locus[locus])}, J={len(j_seqs_by_locus[locus])}")
    d_seqs = read_germline_fasta(os.path.join(BASE_DIR, D_FASTA))
    print(f"  IGH: D={len(d_seqs)}")

    print("Читаю отфильтрованный файл...")
    df = pd.read_csv(INPUT_FILE, sep="\t")
    total = len(df)

    v_dir = build_folder(OUT_ROOT, "v")
    d_dir = build_folder(OUT_ROOT, "d")
    j_dir = build_folder(OUT_ROOT, "j")
    vj_dir = build_folder(OUT_ROOT, "vj")

    v_groups = {}
    d_groups = {}
    j_groups = {}
    vj_groups = {}

    skipped_locus = 0
    for i, (_, row) in enumerate(df.iterrows(), 1):
        seq_id = row["sequence_id"]
        # sequence_vdj — вырезанный участок V(start)..J(end), без лидера и
        # константного региона: сравнивать с V/D/J-гермлайнами нужно именно его,
        # а не полный "сырой" контиг (иначе выравнивание "размазывается").
        seq = row["sequence_vdj"]
        locus = row.get("locus")

        v_seqs = v_seqs_by_locus.get(locus)
        j_seqs = j_seqs_by_locus.get(locus)
        if not isinstance(seq, str) or not seq or v_seqs is None or j_seqs is None:
            skipped_locus += 1
            continue

        best_v, _ = best_germline_match(seq, v_seqs)
        best_j, _ = best_germline_match(seq, j_seqs)

        v_groups.setdefault(best_v, []).append((seq_id, seq))
        j_groups.setdefault(best_j, []).append((seq_id, seq))

        vj_key = f"{best_v}_{best_j}"
        vj_groups.setdefault(vj_key, []).append((seq_id, seq))

        # D-сегмент есть только у тяжёлой цепи (IGH)
        if locus == "IGH":
            best_d, _ = best_germline_match(seq, d_seqs)
            d_groups.setdefault(best_d, []).append((seq_id, seq))

        if i % 50 == 0 or i == total:
            print(f"  обработано {i}/{total}")

    if skipped_locus:
        print(f"Пропущено строк (нет sequence_vdj или неизвестный locus): {skipped_locus}")

    for gene, records in v_groups.items():
        write_fasta(records, v_dir, gene)
    for gene, records in d_groups.items():
        write_fasta(records, d_dir, gene)
    for gene, records in j_groups.items():
        write_fasta(records, j_dir, gene)
    for gene, records in vj_groups.items():
        write_fasta(records, vj_dir, gene)

    print(f"Готово! Fasta-файлов: V={len(v_groups)}, D={len(d_groups)}, J={len(j_groups)}, V+J={len(vj_groups)}")
    print(f"Папки лежат в {OUT_ROOT}/v, /d, /j, /vj")


if __name__ == "__main__":
    main()
