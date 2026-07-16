"""
Тесты для scripts/group_by_germlines.py

Запуск:  pytest tests/test_group_by_germlines.py -v
"""
import os
import sys

import pandas as pd
import pytest

SCRIPTS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts")
sys.path.insert(0, SCRIPTS_DIR)

import group_by_germlines as gg  # noqa: E402


# ---------- read_germline_fasta ----------

def test_read_germline_fasta_basic(tmp_path):
    fasta = tmp_path / "V.fasta"
    fasta.write_text(
        ">acc|IGHV1-2*01|extra\n"
        "acgt\n"
        "acgt\n"
        ">acc2|IGHV1-3*01|extra\n"
        "tttt\n"
    )
    seqs = gg.read_germline_fasta(str(fasta))
    # последовательности приводятся к верхнему регистру, строки склеиваются
    assert seqs == {"IGHV1-2*01": "ACGTACGT", "IGHV1-3*01": "TTTT"}


def test_read_germline_fasta_no_pipe_uses_whole_header(tmp_path):
    fasta = tmp_path / "V.fasta"
    fasta.write_text(">IGHV1-2*01\nacgt\n")
    seqs = gg.read_germline_fasta(str(fasta))
    assert seqs == {"IGHV1-2*01": "ACGT"}


# ---------- safe_filename ----------

@pytest.mark.parametrize("raw,expected", [
    ("IGHV1-2*01", "IGHV1-2_01"),
    ("simple_name.fasta", "simple_name.fasta"),
    ("gene/with/slashes", "gene_with_slashes"),
    ("gene with spaces", "gene_with_spaces"),
])
def test_safe_filename(raw, expected):
    assert gg.safe_filename(raw) == expected


# ---------- best_germline_match ----------

def test_best_germline_match_finds_identical_sequence():
    germline_dict = {
        "GENE_A": "ACGTACGTACGT",
        "GENE_B": "TTTTTTTTTTTT",
    }
    best_name, best_score = gg.best_germline_match("ACGTACGTACGT", germline_dict)
    assert best_name == "GENE_A"
    assert best_score > 0


def test_best_germline_match_picks_closer_sequence():
    germline_dict = {
        "GENE_A": "ACGTACGTACGT",  # похоже
        "GENE_B": "TTTTTTTTTTTT",  # непохоже
    }
    query = "ACGTACGTACGA"  # почти идентично GENE_A, отличается 1 буквой
    best_name, _ = gg.best_germline_match(query, germline_dict)
    assert best_name == "GENE_A"


def test_best_germline_match_case_insensitive():
    germline_dict = {"GENE_A": "ACGTACGT"}
    best_name, best_score = gg.best_germline_match("acgtacgt", germline_dict)
    assert best_name == "GENE_A"
    assert best_score > 0


# ---------- build_folder ----------

def test_build_folder_creates_directory(tmp_path):
    out_root = str(tmp_path / "grouped")
    result_path = gg.build_folder(out_root, "v")
    assert os.path.isdir(result_path)
    assert result_path == os.path.join(out_root, "v")


def test_build_folder_idempotent(tmp_path):
    out_root = str(tmp_path / "grouped")
    gg.build_folder(out_root, "v")
    # повторный вызов не должен падать, даже если папка уже есть
    result_path = gg.build_folder(out_root, "v")
    assert os.path.isdir(result_path)


# ---------- write_fasta ----------

def test_write_fasta_content(tmp_path):
    records = [("seq1", "ACGT"), ("seq2", "TTTT")]
    gg.write_fasta(records, str(tmp_path), "IGHV1-2*01")

    fpath = tmp_path / "IGHV1-2_01.fasta"
    assert fpath.exists()
    content = fpath.read_text()
    assert content == ">seq1\nACGT\n>seq2\nTTTT\n"


def test_write_fasta_uses_safe_filename(tmp_path):
    gg.write_fasta([("seq1", "ACGT")], str(tmp_path), "gene/with slashes")
    fpath = tmp_path / "gene_with_slashes.fasta"
    assert fpath.exists()


# ---------- end-to-end через main() на маленьком наборе данных ----------

def test_main_end_to_end(tmp_path, monkeypatch):
    # маленькие germline-справочники для всех локусов + D для IGH
    for locus, files in gg.LOCUS_VJ_FASTA.items():
        (tmp_path / files["v"]).write_text(">acc|V1*01|x\nACGTACGTACGT\n")
        (tmp_path / files["j"]).write_text(">acc|J1*01|x\nTTTTTTTT\n")
    (tmp_path / gg.D_FASTA).write_text(">acc|D1*01|x\nGGGG\n")

    filtered_tsv = tmp_path / "BCR_data_filtered.tsv"
    seq = "ACGTACGTACGTTTTTTTTT"
    filtered_tsv.write_text(
        "sequence_id\tsequence_vdj\tlocus\n"
        f"read1\t{seq}\tIGH\n"
        f"read2\t{seq}\tIGK\n"
    )

    out_root = tmp_path / "grouped_by_germlines"
    monkeypatch.setattr(gg, "BASE_DIR", str(tmp_path))
    monkeypatch.setattr(gg, "INPUT_FILE", str(filtered_tsv))
    monkeypatch.setattr(gg, "OUT_ROOT", str(out_root))

    gg.main()

    assert (out_root / "v").is_dir()
    assert (out_root / "j").is_dir()
    assert (out_root / "vj").is_dir()
    assert (out_root / "d").is_dir()

    # у IGH-строки должен появиться файл в папке d (D-сегмент только у тяжёлой цепи)
    d_files = list((out_root / "d").iterdir())
    assert len(d_files) == 1

    # в v/j/vj должно быть суммарно 2 обработанные записи (IGH + IGK)
    v_files = list((out_root / "v").iterdir())
    assert len(v_files) >= 1


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
