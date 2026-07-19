"""
Тесты для scripts/filter_by_symbol_count/filter_by_symbol_count.py

Запуск:  pytest tests/test_filter_by_symbol_count.py -v
"""
import os
import sys

import pytest

SCRIPTS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts")
sys.path.insert(0, SCRIPTS_DIR)

import filter_by_symbol_count as fbc  # noqa: E402


# ---------- count_symbol ----------

def test_count_symbol_basic(tmp_path):
    fasta = tmp_path / "genes.fasta"
    fasta.write_text(">seq1\nACGT\n>seq2\nTTTT\n>seq3\nGGGG\n")
    assert fbc.count_symbol(str(fasta)) == 3


def test_count_symbol_no_matches(tmp_path):
    fasta = tmp_path / "empty.fasta"
    fasta.write_text("ACGTACGT\n")
    assert fbc.count_symbol(str(fasta)) == 0


def test_count_symbol_custom_symbol(tmp_path):
    f = tmp_path / "data.txt"
    f.write_text("a#b#c")
    assert fbc.count_symbol(str(f), symbol="#") == 2


# ---------- filter_and_copy_files ----------

def _make_fasta(path, n_sequences):
    with open(path, "w", encoding="utf-8") as f:
        for i in range(n_sequences):
            f.write(f">seq{i}\nACGT\n")


def test_filter_and_copy_files_keeps_only_in_range(tmp_path):
    target_dir = tmp_path / "vj"
    output_dir = tmp_path / "vj_filtered"
    target_dir.mkdir()

    _make_fasta(target_dir / "too_small.fasta", 2)     # <=5 -> не проходит
    _make_fasta(target_dir / "good_10.fasta", 10)       # проходит
    _make_fasta(target_dir / "good_50.fasta", 50)        # проходит
    _make_fasta(target_dir / "too_big.fasta", 150)       # >=100 -> не проходит

    copied = fbc.filter_and_copy_files(str(target_dir), str(output_dir))

    assert copied == 2
    result_files = sorted(os.listdir(output_dir))
    assert result_files == ["good_10.fasta", "good_50.fasta"]


def test_filter_and_copy_files_boundary_values_excluded(tmp_path):
    target_dir = tmp_path / "vj"
    output_dir = tmp_path / "vj_filtered"
    target_dir.mkdir()

    _make_fasta(target_dir / "exactly_5.fasta", 5)     # граница не включена
    _make_fasta(target_dir / "exactly_100.fasta", 100)  # граница не включена

    copied = fbc.filter_and_copy_files(str(target_dir), str(output_dir))

    assert copied == 0


def test_filter_and_copy_files_creates_output_dir(tmp_path):
    target_dir = tmp_path / "vj"
    output_dir = tmp_path / "new_folder" / "vj_filtered"
    target_dir.mkdir()
    _make_fasta(target_dir / "good.fasta", 10)

    fbc.filter_and_copy_files(str(target_dir), str(output_dir))

    assert output_dir.is_dir()


def test_filter_and_copy_files_empty_folder(tmp_path):
    target_dir = tmp_path / "vj"
    output_dir = tmp_path / "vj_filtered"
    target_dir.mkdir()

    copied = fbc.filter_and_copy_files(str(target_dir), str(output_dir))

    assert copied == 0
    assert output_dir.is_dir()


# ---------- end-to-end через main() ----------
#
# Как и в остальных скриптах, main() использует get_paths() из paths.py,
# поэтому подменяем саму функцию get_paths, а не жёстко прописанные пути.

def test_main_end_to_end(tmp_path, monkeypatch):
    output_dir = tmp_path / "results" / "testkey"
    target_dir = output_dir / "grouped_by_germlines" / "vj"
    target_dir.mkdir(parents=True)

    _make_fasta(target_dir / "too_small.fasta", 2)
    _make_fasta(target_dir / "good_10.fasta", 10)
    _make_fasta(target_dir / "too_big.fasta", 150)

    def fake_get_paths(key, create_output=True):
        return {
            "input_dir": str(tmp_path / "data" / key),
            "output_dir": str(output_dir),
            "imgt_dir": str(tmp_path / "data" / "IMGT" / "Homo_sapiens" / "IG"),
        }

    monkeypatch.setattr(fbc, "get_paths", fake_get_paths)
    monkeypatch.setattr(sys, "argv", ["filter_by_symbol_count.py", "--key", "testkey"])

    fbc.main()

    # новая папка теперь на уровне results/<key>/, а НЕ рядом с vj внутри grouped_by_germlines
    filtered_dir = output_dir / "vj_filtered"
    assert filtered_dir.is_dir()
    result_files = sorted(os.listdir(filtered_dir))
    assert result_files == ["good_10.fasta"]


def test_main_end_to_end_missing_target_dir_does_not_crash(tmp_path, monkeypatch, capsys):
    output_dir = tmp_path / "results" / "testkey"
    output_dir.mkdir(parents=True)  # но подпапки grouped_by_germlines/vj нет

    def fake_get_paths(key, create_output=True):
        return {
            "input_dir": str(tmp_path / "data" / key),
            "output_dir": str(output_dir),
            "imgt_dir": str(tmp_path / "data" / "IMGT" / "Homo_sapiens" / "IG"),
        }

    monkeypatch.setattr(fbc, "get_paths", fake_get_paths)
    monkeypatch.setattr(sys, "argv", ["filter_by_symbol_count.py", "--key", "testkey"])

    fbc.main()  # не должно падать с исключением

    captured = capsys.readouterr()
    assert "не найден" in captured.out


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
