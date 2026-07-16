import argparse
import os
import sys
from collections import Counter

import pandas as pd

# ============ ПУТИ ============
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# paths.py лежит на уровень выше (в scripts/), а не рядом с этим файлом,
# поэтому добавляем scripts/ в sys.path, чтобы его можно было импортировать
SCRIPTS_DIR = os.path.dirname(BASE_DIR)
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from paths import get_paths  # noqa: E402

MIN_JUNCTION = 15
MAX_JUNCTION = 100
V_MIN_FRACTION = 0.5


def read_germline_lengths(path):
    lengths = {}
    name = None
    chunks = []

    def flush():
        if name is not None:
            lengths[name] = len("".join(chunks))

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                parts = line[1:].split("|")
                name = parts[1].strip() if len(parts) > 1 else parts[0].strip()
                chunks = []
            else:
                chunks.append(line)
        flush()
    return lengths


def primary_gene(call_str):
    """Берёт первый ген из списка (через запятую) и убирает суффиксы вида
    '_T194C' (точечные мутации, которые некоторые пайплайны дописывают к
    названию аллели) — без этого имя не совпадает со справочником."""
    if pd.isna(call_str) or not str(call_str).strip():
        return None
    gene = str(call_str).split(",")[0].strip()
    gene = gene.split("_")[0]
    return gene or None


def lookup_length(gene_name, lengths_dict):
    """Ищет длину гена. Сначала точное совпадение, потом — без учёта аллели (*01, *02...)."""
    if gene_name is None:
        return None
    if gene_name in lengths_dict:
        return lengths_dict[gene_name]

    # fallback: искать без учёта конкретной аллели (до символа *)
    base = gene_name.split("*")[0]
    candidates = [v for k, v in lengths_dict.items() if k.split("*")[0] == base]
    if candidates:
        return sum(candidates) / len(candidates)  # средняя длина по всем аллелям гена
    return None


def read_input(path):
    ext = os.path.splitext(path)[1].lower()
    sep = "\t" if ext == ".tsv" else ","
    return pd.read_csv(path, sep=sep)


def check_row(row, v_lengths_by_locus, j_lengths_by_locus, reasons):
    locus = row.get("locus")
    v_lengths = v_lengths_by_locus.get(locus)
    j_lengths = j_lengths_by_locus.get(locus)
    if v_lengths is None or j_lengths is None:
        reasons["неизвестный/неподдерживаемый locus"] += 1
        return False

    v_gene = primary_gene(row.get("v_call"))
    j_gene = primary_gene(row.get("j_call"))
    # sequence_vdj — это уже вырезанный участок V(start)..J(end), без лидерной
    # последовательности и константного региона. Полная "sequence" (raw contig)
    # всегда длиннее v_len+j_len+junction, поэтому по ней проверка длины не
    # проходит почти никогда — сравнивать нужно именно с sequence_vdj.
    seq = row.get("sequence_vdj")

    if not isinstance(seq, str) or not seq:
        reasons["нет последовательности"] += 1
        return False
    if v_gene is None or j_gene is None:
        reasons["нет v_call/j_call"] += 1
        return False

    v_len = lookup_length(v_gene, v_lengths)
    j_len = lookup_length(j_gene, j_lengths)
    if v_len is None:
        reasons["V-ген не найден в справочнике"] += 1
        return False
    if j_len is None:
        reasons["J-ген не найден в справочнике"] += 1
        return False

    min_expected = V_MIN_FRACTION * v_len + MIN_JUNCTION + j_len
    max_expected = v_len + MAX_JUNCTION + j_len

    if not (min_expected <= len(seq) <= max_expected):
        reasons["длина вне ожидаемого диапазона"] += 1
        return False

    return True


def parse_args():
    parser = argparse.ArgumentParser(description="Фильтрация BCR-последовательностей.")
    parser.add_argument(
        "-k", "--key",
        required=True,
        help="Название подпапки внутри data/ и results/ (например: batch1)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    paths = get_paths(args.key)

    input_dir = paths["input_dir"]
    output_dir = paths["output_dir"]
    imgt_dir = paths["imgt_dir"]  # data/IMGT/Homo_sapiens/IG — подаётся автоматически

    input_file = os.path.join(input_dir, "BCR_data.tsv")
    output_file = os.path.join(output_dir, "BCR_data_filtered.tsv")

    # germline-справочники лежат в imgt_dir (см. paths.py)
    locus_fasta = {
        "IGH": {"v": os.path.join(imgt_dir, "IGHV.fasta"), "j": os.path.join(imgt_dir, "IGHJ.fasta")},
        "IGK": {"v": os.path.join(imgt_dir, "IGKV.fasta"), "j": os.path.join(imgt_dir, "IGKJ.fasta")},
        "IGL": {"v": os.path.join(imgt_dir, "IGLV.fasta"), "j": os.path.join(imgt_dir, "IGLJ.fasta")},
    }

    print("Читаю germline-справочники (только длины)...")
    v_lengths_by_locus = {}
    j_lengths_by_locus = {}
    for locus, files in locus_fasta.items():
        v_lengths_by_locus[locus] = read_germline_lengths(files["v"])
        j_lengths_by_locus[locus] = read_germline_lengths(files["j"])
        print(f"  {locus}: V={len(v_lengths_by_locus[locus])}, J={len(j_lengths_by_locus[locus])}")

    print(f"Читаю входной файл {input_file} ...")
    df = read_input(input_file)
    total_before = len(df)
    print("Пример значений v_call в твоих данных:", df["v_call"].dropna().unique()[:5])
    print("Локусы в данных:", df["locus"].value_counts().to_dict())

    df = df.drop_duplicates(subset="sequence", keep="first")
    after_dedup = len(df)

    reasons = Counter()
    mask = df.apply(lambda row: check_row(row, v_lengths_by_locus, j_lengths_by_locus, reasons), axis=1)
    df_valid = df[mask]
    after_length_filter = len(df_valid)

    print(f"\nВсего было: {total_before}")
    print(f"После удаления дубликатов: {after_dedup} (убрано {total_before - after_dedup})")
    print(f"После проверки длины: {after_length_filter} (убрано {after_dedup - after_length_filter})")
    print("\nПричины отбраковки:")
    for reason, count in reasons.most_common():
        print(f"  {reason}: {count}")

    df_valid.to_csv(output_file, sep="\t", index=False)
    print(f"\nГотово! Отфильтрованные данные сохранены в {output_file}")


if __name__ == "__main__":
    main()
