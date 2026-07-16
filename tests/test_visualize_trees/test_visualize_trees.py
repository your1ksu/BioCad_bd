import subprocess
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]
SCRIPT = BASE / "scripts" / "visualize_trees" / "visualize_trees.sh"
DATA_DIR = Path(__file__).parent / "test_visualize_trees_data"


def test_successful_run(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    result_dir = tmp_path / "test_result"

    result = subprocess.run(
        ["bash", str(SCRIPT), str(DATA_DIR), str(result_dir)],
        capture_output=True, text=True, timeout=300,
    )

    assert result.returncode == 0
    assert result_dir.is_dir()

    for family_dir in DATA_DIR.iterdir():
        if not family_dir.is_dir():
            continue
        name = family_dir.name
        subdir = result_dir / name
        assert subdir.is_dir(), f"{subdir} not created"
        assert (subdir / f"{name}.html").exists(), f"html for {name} missing"


def test_nonexistent_input_dir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fake_dir = tmp_path / "does_not_exist"
    result_dir = tmp_path / "test_result"

    result = subprocess.run(
        ["bash", str(SCRIPT), str(fake_dir), str(result_dir)],
        capture_output=True, text=True,
    )

    assert result.returncode != 0
    assert "not found" in result.stderr or "not found" in result.stdout