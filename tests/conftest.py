import os
import subprocess
import sys

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(TESTS_DIR)
SCRIPTS_DIR = os.path.join(REPO_ROOT, "scripts")

# добавляем саму scripts/, а также все её подпапки первого уровня —
# на случай, если каждый скрипт лежит в своей отдельной подпапке
# (например scripts/filtered/, scripts/group_by_germlines/)
paths_to_add = [SCRIPTS_DIR]
if os.path.isdir(SCRIPTS_DIR):
    for name in os.listdir(SCRIPTS_DIR):
        full_path = os.path.join(SCRIPTS_DIR, name)
        if os.path.isdir(full_path):
            paths_to_add.append(full_path)

for path in paths_to_add:
    if path not in sys.path:
        sys.path.insert(0, path)

ENV_NAME = "biocad_bcr_pipeline_environment"


def pipeline_python():
    """Return a command prefix to run Python inside the pipeline conda environment."""
    try:
        subprocess.run(
            ["conda", "run", "-n", ENV_NAME, "python", "-c", "import Bio"],
            capture_output=True, text=True, check=True
        )
        return ["conda", "run", "-n", ENV_NAME, "python"]
    except Exception:
        return [sys.executable]
