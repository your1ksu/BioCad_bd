import subprocess
from pathlib import Path

import conftest

BASE = Path(__file__).resolve().parents[2]
SCRIPT = BASE / "scripts" / "04a_build_trees_iqtree" / "build_trees_iqtree.py"
DATA_DIR = Path(__file__).parent / "test_build_iqtree_data"


def test_successful_run(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    result_dir = tmp_path / "test_result"

    cmd = conftest.pipeline_python() + [str(SCRIPT), "-i", str(DATA_DIR), "-o", str(result_dir)]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

    assert result.returncode == 0, result.stderr
    assert result_dir.is_dir()

    for fasta in DATA_DIR.glob("*.fasta"):
        name = fasta.stem
        subdir = result_dir / name
        assert subdir.is_dir(), f"{subdir} not created"
        assert (subdir / f"{name}.treefile").exists(), f"treefile for {name} missing"


def test_nonexistent_input_dir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fake_dir = tmp_path / "does_not_exist"
    result_dir = tmp_path / "test_result"

    cmd = conftest.pipeline_python() + [str(SCRIPT), "-i", str(fake_dir), "-o", str(result_dir)]
    result = subprocess.run(cmd, capture_output=True, text=True)

    assert result.returncode != 0
