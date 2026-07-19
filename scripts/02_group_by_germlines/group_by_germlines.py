#!/usr/bin/env python3
"""
Группировка BCR-последовательностей по germline-генам.

Использование:
    python group_by_germlines.py \
        -i results/report_17072026_143000/BCR_data_filtered.tsv \
        -o results/report_17072026_143000/grouped_by_germlines \
        -r data/references
"""

import argparse
import os
import sys
from pathlib import Path

import pandas as pd
from Bio import Align

_aligner = Align.PairwiseAligner()
_aligner.mode = "local"

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from shared.utils import read_germline_fasta, safe_filename, format_group_key





def parse_vj_call(call_str, germline_dict):
    """Парсит v_call/j_call из AIRR (например 'IGHV1-2*01') и ищет в germline_dict.
    Возвращает точное совпадение или None."""
    if pd.isna(call_str) or not str(call_str).strip():
        return None
    # Берём первый ген до запятой
    gene = str(call_str).split(",")[0].strip()
    gene = gene.split("_")[0]  # убираем суффиксы типа _T194C
    if gene in germline_dict:
        return gene
    # Пробуем без аллели (*01, *02...)
    base = gene.split("*")[0]
    for k in germline_dict:
        if k.startswith(base + "*") or k == base:
            return k
    return None


def best_germline_match(query_seq, germline_dict):
    """Находит ген с максимальным скором выравнивания (наиболее похожий гермлайн)."""
    best_name = None
    best_score = float("-inf")
    query_seq = query_seq.upper()
    for name, gseq in germline_dict.items():
        score = _aligner.score(query_seq, gseq)
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


def parse_args():
    parser = argparse.ArgumentParser(description="Группировка BCR-последовательностей по germline-генам.")
    parser.add_argument("-i", "--input", required=True,
                        help="Путь к входному TSV файлу (результат filter_sequences.py).")
    parser.add_argument("-o", "--output", required=True,
                        help="Путь к выходной директории (будет создана grouped_by_germlines внутри).")
    parser.add_argument("-r", "--ref-dir", required=True,
                        help="Папка с гермлайновыми справочниками.")
    parser.add_argument("--grouping-strategy", choices=["allele", "gene", "v_only"], default="gene",
                        help="Стратегия группировки (default: gene)")
    return parser.parse_args()


def main():
    args = parse_args()
    input_file = args.input
    output_dir = args.output
    ref_dir = args.ref_dir
    strategy = args.grouping_strategy

    if not os.path.isfile(input_file):
        print(f"Ошибка: входной файл не найден: {input_file}", file=sys.stderr)
        sys.exit(1)
    if not os.path.isdir(ref_dir):
        print(f"Ошибка: папка референсов не найдена: {ref_dir}", file=sys.stderr)
        sys.exit(1)

    out_root = output_dir

    locus_vj_fasta = {
        "IGH": {"v": os.path.join(ref_dir, "IGHV.fasta"), "j": os.path.join(ref_dir, "IGHJ.fasta")},
        "IGK": {"v": os.path.join(ref_dir, "IGKV.fasta"), "j": os.path.join(ref_dir, "IGKJ.fasta")},
        "IGL": {"v": os.path.join(ref_dir, "IGLV.fasta"), "j": os.path.join(ref_dir, "IGLJ.fasta")},
    }

    print("Читаю germline-справочники...")
    v_seqs_by_locus = {}
    j_seqs_by_locus = {}
    for locus, files in locus_vj_fasta.items():
        v_seqs_by_locus[locus] = read_germline_fasta(files["v"])
        j_seqs_by_locus[locus] = read_germline_fasta(files["j"])
        print(f"  {locus}: V={len(v_seqs_by_locus[locus])}, J={len(j_seqs_by_locus[locus])}")

    print(f"Читаю отфильтрованный файл {input_file} ...")
    df = pd.read_csv(input_file, sep="\t")
    total = len(df)

    vj_dir = build_folder(out_root, "vj")
    vj_groups = {}

    skipped_locus = 0
    aligned_v = 0
    aligned_j = 0

    for i, (_, row) in enumerate(df.iterrows(), 1):
        seq_id = row["sequence_id"]
        seq = row["sequence_vdj"]
        locus = row.get("locus")

        v_seqs = v_seqs_by_locus.get(locus)
        j_seqs = j_seqs_by_locus.get(locus)
        if not isinstance(seq, str) or not seq or v_seqs is None or j_seqs is None:
            skipped_locus += 1
            continue

        v_call = row.get("v_call")
        j_call = row.get("j_call")

        best_v = parse_vj_call(v_call, v_seqs)
        if best_v is None:
            best_v, _ = best_germline_match(seq, v_seqs)
            aligned_v += 1

        best_j = parse_vj_call(j_call, j_seqs)
        if best_j is None:
            best_j, _ = best_germline_match(seq, j_seqs)
            aligned_j += 1

        vj_key = format_group_key(best_v, best_j, strategy)
        vj_groups.setdefault(vj_key, []).append((seq_id, seq))

        if i % 500 == 0 or i == total:
            print(f"  обработано {i}/{total}")

    if aligned_v or aligned_j:
        print(f"  Выравнивание использовано для {aligned_v} V-генов и {aligned_j} J-генов (остальные по v_call/j_call)")
    if skipped_locus:
        print(f"Пропущено строк (нет sequence_vdj или неизвестный locus): {skipped_locus}")

    for gene, records in vj_groups.items():
        write_fasta(records, vj_dir, gene)

    print(f"Готово! Fasta-файлов: V+J={len(vj_groups)}")
    print(f"Стратегия группировки: {strategy}")
    print(f"Папка: {out_root}/vj")


if __name__ == "__main__":
    main()