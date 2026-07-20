import os
import sys
from pathlib import Path

import pytest

SCRIPTS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts")
sys.path.insert(0, SCRIPTS_DIR)

import importlib.util
spec = importlib.util.spec_from_file_location(
    "filter_by_symbol_count",
    os.path.join(SCRIPTS_DIR, "filter_by_symbol_count", "filter_by_symbol_count.py")
)
fbc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fbc)


def test_count_sequences_basic(tmp_path):
    fasta = tmp_path / "genes.fasta"
    fasta.write_text(">seq1\nACGT\n>seq2\nTTTT\n>seq3\nGGGG\n")
    assert fbc.count_sequences(fasta) == 3


def test_count_sequences_no_matches(tmp_path):
    fasta = tmp_path / "empty.fasta"
    fasta.write_text("ACGTACGT\n")
    assert fbc.count_sequences(fasta) == 0


def _make_fasta(path, n_sequences):
    with open(path, "w", encoding="utf-8") as f:
        for i in range(n_sequences):
            f.write(f">seq{i}\nACGT\n")


def test_filter_and_link_files_keeps_only_in_range(tmp_path):
    target_dir = tmp_path / "vj"
    output_dir = tmp_path / "vj_filtered"
    target_dir.mkdir()

    _make_fasta(target_dir / "too_small.fasta", 2)
    _make_fasta(target_dir / "good_10.fasta", 10)
    _make_fasta(target_dir / "good_50.fasta", 50)
    _make_fasta(target_dir / "too_big.fasta", 150)

    copied = fbc.filter_and_link_files(target_dir, output_dir, min_count=5, max_count=100)

    assert copied == 2
    result_files = sorted(os.listdir(output_dir))
    assert result_files == ["good_10.fasta", "good_50.fasta"]


def test_filter_and_link_files_boundary_values_excluded(tmp_path):
    target_dir = tmp_path / "vj"
    output_dir = tmp_path / "vj_filtered"
    target_dir.mkdir()

    _make_fasta(target_dir / "exactly_5.fasta", 5)
    _make_fasta(target_dir / "exactly_100.fasta", 100)

    copied = fbc.filter_and_link_files(target_dir, output_dir, min_count=5, max_count=100)

    assert copied == 0


def test_filter_and_link_files_creates_output_dir(tmp_path):
    target_dir = tmp_path / "vj"
    output_dir = tmp_path / "new_folder" / "vj_filtered"
    target_dir.mkdir()
    _make_fasta(target_dir / "good.fasta", 10)

    fbc.filter_and_link_files(target_dir, output_dir, min_count=5, max_count=100)

    assert output_dir.is_dir()


def test_filter_and_link_files_empty_folder(tmp_path):
    target_dir = tmp_path / "vj"
    output_dir = tmp_path / "vj_filtered"
    target_dir.mkdir()

    copied = fbc.filter_and_link_files(target_dir, output_dir, min_count=5, max_count=100)

    assert copied == 0
    assert output_dir.is_dir()


def test_main_end_to_end(tmp_path, monkeypatch):
    target_dir = tmp_path / "vj"
    target_dir.mkdir(parents=True)

    _make_fasta(target_dir / "too_small.fasta", 2)
    _make_fasta(target_dir / "good_10.fasta", 10)
    _make_fasta(target_dir / "too_big.fasta", 150)

    output_dir = tmp_path / "vj_filtered"
    monkeypatch.setattr(sys, "argv", [
        "filter_by_symbol_count.py",
        "-i", str(target_dir),
        "-o", str(output_dir),
        "--min", "5",
        "--max", "100",
    ])

    fbc.main()

    assert output_dir.is_dir()
    result_files = sorted(os.listdir(output_dir))
    assert result_files == ["good_10.fasta"]


def test_main_end_to_end_missing_target_dir_does_not_crash(tmp_path, monkeypatch, capsys):
    target_dir = tmp_path / "nonexistent"
    output_dir = tmp_path / "vj_filtered"

    monkeypatch.setattr(sys, "argv", [
        "filter_by_symbol_count.py",
        "-i", str(target_dir),
        "-o", str(output_dir),
    ])

    with pytest.raises(SystemExit):
        fbc.main()

    captured = capsys.readouterr()
    assert "не найдена" in captured.out or "не найдена" in captured.err


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
