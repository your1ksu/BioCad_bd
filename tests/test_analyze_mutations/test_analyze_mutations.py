import subprocess
from pathlib import Path

import conftest

BASE = Path(__file__).resolve().parents[2]
SCRIPT = BASE / "scripts" / "06_analyze_mutations" / "run_mutations.py"
DATA_DIR = Path(__file__).parent / "test_analyze_mutations_data"
REF_DIR = BASE / "data" / "references"


def test_successful_run(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    result_dir = tmp_path / "test_result"

    cmd = conftest.pipeline_python() + [str(SCRIPT), "-i", str(DATA_DIR), "-o", str(result_dir), "-r", str(REF_DIR)]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

    assert result.returncode == 0, result.stderr
    assert result_dir.is_dir()

    for fasta in DATA_DIR.glob("*.fasta"):
        name = fasta.stem
        subdir = result_dir / name
        assert subdir.is_dir(), f"{subdir} not created"
        assert (subdir / "mutations.tsv").exists(), f"mutations.tsv for {name} missing"
        assert (subdir / "mutations_summary.tsv").exists(), f"mutations_summary.tsv for {name} missing"


def test_nonexistent_input_dir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fake_dir = tmp_path / "does_not_exist"
    result_dir = tmp_path / "test_result"

    cmd = conftest.pipeline_python() + [str(SCRIPT), "-i", str(fake_dir), "-o", str(result_dir), "-r", str(REF_DIR)]
    result = subprocess.run(cmd, capture_output=True, text=True)

    assert result.returncode != 0
