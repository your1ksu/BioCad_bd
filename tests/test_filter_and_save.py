"""
Тесты для scripts/filter_and_save.py

Запуск:  pytest tests/test_filter_and_save.py -v
"""
import os
import sys
from collections import Counter

import pandas as pd
import pytest

# scripts/ не является пакетом (нет __init__.py), поэтому добавляем путь вручную
SCRIPTS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts")
sys.path.insert(0, SCRIPTS_DIR)

import filter_and_save as fs  # noqa: E402


# ---------- primary_gene ----------

def test_primary_gene_simple():
    assert fs.primary_gene("IGHV1-2*01") == "IGHV1-2*01"


def test_primary_gene_takes_first_of_several():
    # некоторые пайплайны пишут несколько генов через запятую — берём первый
    assert fs.primary_gene("IGHV1-2*01,IGHV1-2*02") == "IGHV1-2*01"


def test_primary_gene_strips_mutation_suffix():
    # суффикс вида _T194C (точечная мутация) не входит в имя гена
    assert fs.primary_gene("IGHV1-2*01_T194C") == "IGHV1-2*01"


def test_primary_gene_none_for_empty_or_nan():
    assert fs.primary_gene(None) is None
    assert fs.primary_gene("") is None
    assert fs.primary_gene(float("nan")) is None


# ---------- read_germline_lengths ----------

def test_read_germline_lengths_basic(tmp_path):
    fasta = tmp_path / "V.fasta"
    fasta.write_text(
        ">accession|IGHV1-2*01|extra\n"
        "ACGTACGT\n"
        "ACGT\n"
        ">accession2|IGHV1-2*02|extra\n"
        "ACGTACGTACGT\n"
    )
    lengths = fs.read_germline_lengths(str(fasta))
    assert lengths == {"IGHV1-2*01": 12, "IGHV1-2*02": 12}


def test_read_germline_lengths_no_pipe_uses_whole_header(tmp_path):
    fasta = tmp_path / "V.fasta"
    fasta.write_text(">IGHV1-2*01\nACGT\n")
    lengths = fs.read_germline_lengths(str(fasta))
    assert lengths == {"IGHV1-2*01": 4}


# ---------- lookup_length ----------

def test_lookup_length_exact_match():
    lengths = {"IGHV1-2*01": 100}
    assert fs.lookup_length("IGHV1-2*01", lengths) == 100


def test_lookup_length_fallback_averages_alleles():
    lengths = {"IGHV1-2*01": 100, "IGHV1-2*02": 102}
    # запрошенной аллели *03 нет -> усредняем по всем аллелям гена IGHV1-2
    assert fs.lookup_length("IGHV1-2*03", lengths) == 101


def test_lookup_length_none_when_gene_not_found():
    lengths = {"IGHV1-2*01": 100}
    assert fs.lookup_length("IGHV9-9*01", lengths) is None


def test_lookup_length_none_for_none_gene():
    assert fs.lookup_length(None, {"IGHV1-2*01": 100}) is None


# ---------- check_row ----------

def _make_lookups():
    v_lengths_by_locus = {"IGH": {"IGHV1-2*01": 100}}
    j_lengths_by_locus = {"IGH": {"IGHJ4*01": 20}}
    return v_lengths_by_locus, j_lengths_by_locus


def test_check_row_valid_passes():
    v_lengths, j_lengths = _make_lookups()
    # ожидаемая длина: 0.5*100 + 15 .. 100 + 100 + 20 = [65, 220]
    seq = "A" * 150
    row = pd.Series({
        "locus": "IGH",
        "v_call": "IGHV1-2*01",
        "j_call": "IGHJ4*01",
        "sequence_vdj": seq,
    })
    reasons = Counter()
    assert fs.check_row(row, v_lengths, j_lengths, reasons) is True
    assert len(reasons) == 0


def test_check_row_unknown_locus_rejected():
    v_lengths, j_lengths = _make_lookups()
    row = pd.Series({
        "locus": "TCR",  # не поддерживается
        "v_call": "IGHV1-2*01",
        "j_call": "IGHJ4*01",
        "sequence_vdj": "A" * 150,
    })
    reasons = Counter()
    assert fs.check_row(row, v_lengths, j_lengths, reasons) is False
    assert reasons["неизвестный/неподдерживаемый locus"] == 1


def test_check_row_missing_sequence_rejected():
    v_lengths, j_lengths = _make_lookups()
    row = pd.Series({
        "locus": "IGH",
        "v_call": "IGHV1-2*01",
        "j_call": "IGHJ4*01",
        "sequence_vdj": None,
    })
    reasons = Counter()
    assert fs.check_row(row, v_lengths, j_lengths, reasons) is False
    assert reasons["нет последовательности"] == 1


def test_check_row_missing_v_or_j_call_rejected():
    v_lengths, j_lengths = _make_lookups()
    row = pd.Series({
        "locus": "IGH",
        "v_call": None,
        "j_call": "IGHJ4*01",
        "sequence_vdj": "A" * 150,
    })
    reasons = Counter()
    assert fs.check_row(row, v_lengths, j_lengths, reasons) is False
    assert reasons["нет v_call/j_call"] == 1


def test_check_row_gene_not_in_reference_rejected():
    v_lengths, j_lengths = _make_lookups()
    row = pd.Series({
        "locus": "IGH",
        "v_call": "IGHV9-9*01",  # нет в справочнике
        "j_call": "IGHJ4*01",
        "sequence_vdj": "A" * 150,
    })
    reasons = Counter()
    assert fs.check_row(row, v_lengths, j_lengths, reasons) is False
    assert reasons["V-ген не найден в справочнике"] == 1


def test_check_row_length_out_of_range_rejected():
    v_lengths, j_lengths = _make_lookups()
    row = pd.Series({
        "locus": "IGH",
        "v_call": "IGHV1-2*01",
        "j_call": "IGHJ4*01",
        "sequence_vdj": "A" * 10,  # короче минимума (65)
    })
    reasons = Counter()
    assert fs.check_row(row, v_lengths, j_lengths, reasons) is False
    assert reasons["длина вне ожидаемого диапазона"] == 1


# ---------- read_input ----------

def test_read_input_tsv(tmp_path):
    tsv = tmp_path / "data.tsv"
    tsv.write_text("a\tb\n1\t2\n")
    df = fs.read_input(str(tsv))
    assert list(df.columns) == ["a", "b"]
    assert df.iloc[0]["a"] == 1


def test_read_input_csv(tmp_path):
    csv = tmp_path / "data.csv"
    csv.write_text("a,b\n1,2\n")
    df = fs.read_input(str(csv))
    assert list(df.columns) == ["a", "b"]
    assert df.iloc[0]["b"] == 2


# ---------- end-to-end через main() на маленьком наборе данных ----------

def test_main_end_to_end(tmp_path, monkeypatch):
    # готовим mini germline-справочники для всех локусов, которые ожидает скрипт
    for locus, files in fs.LOCUS_FASTA.items():
        (tmp_path / files["v"]).write_text(">acc|V1*01|extra\n" + "A" * 100 + "\n")
        (tmp_path / files["j"]).write_text(">acc|J1*01|extra\n" + "A" * 20 + "\n")

    # входной файл: 1 валидная строка, 1 дубликат, 1 с неизвестным locus
    input_tsv = tmp_path / "BCR_data.tsv"
    valid_seq = "A" * 150
    input_tsv.write_text(
        "sequence\tsequence_vdj\tlocus\tv_call\tj_call\n"
        f"{valid_seq}\t{valid_seq}\tIGH\tV1*01\tJ1*01\n"
        f"{valid_seq}\t{valid_seq}\tIGH\tV1*01\tJ1*01\n"  # дубликат по sequence
        f"XXXX\tXXXX\tUNKNOWN\tV1*01\tJ1*01\n"
    )

    monkeypatch.setattr(fs, "BASE_DIR", str(tmp_path))
    monkeypatch.setattr(fs, "INPUT_FILE", str(input_tsv))
    output_tsv = tmp_path / "BCR_data_filtered.tsv"
    monkeypatch.setattr(fs, "OUTPUT_FILE", str(output_tsv))

    fs.main()

    assert output_tsv.exists()
    result = pd.read_csv(output_tsv, sep="\t")
    assert len(result) == 1
    assert result.iloc[0]["locus"] == "IGH"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
