#!/usr/bin/env bash

set -euo pipefail

INPUT_DIR="${1:-trees}"
OUTPUT_DIR="${2:-${INPUT_DIR}_visualization}"
CONDA_ENV="trees_building_env"

if ! command -v conda &>/dev/null; then
    echo "Error: conda not found. Install miniconda first."
    exit 1
fi

source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true

if conda env list | grep -q "^${CONDA_ENV} "; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Conda env '${CONDA_ENV}' exists, activating ..."
else
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Creating conda env '${CONDA_ENV}' ..."
    conda create -n "${CONDA_ENV}" -y
fi

conda activate "${CONDA_ENV}"

if ! python -c "import toytree" 2>/dev/null; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Installing toytree in '${CONDA_ENV}' ..."
    conda install -c conda-forge toytree -y
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: directory $INPUT_DIR not found"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python "$script_dir/visualize_trees.py" "$INPUT_DIR" "$OUTPUT_DIR"
