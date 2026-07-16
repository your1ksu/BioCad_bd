#!/usr/bin/env bash

set -euo pipefail

INPUT_DIR="${1:-aligned_sequences}"
OUTPUT_DIR="${2:-trees}"
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

if ! command -v iqtree &>/dev/null; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Installing iqtree in '${CONDA_ENV}' ..."
    conda install -c bioconda iqtree -y
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: directory $INPUT_DIR not found"
    exit 1
fi

if command -v nproc &>/dev/null; then
    NPROC=$(nproc)
elif sysctl -n hw.ncpu &>/dev/null 2>&1; then
    NPROC=$(sysctl -n hw.ncpu)
else
    NPROC=1
fi

mkdir -p "$OUTPUT_DIR"

fasta_files=()
while IFS= read -r -d '' f; do
    fasta_files+=("$f")
done < <(find "$INPUT_DIR" -maxdepth 1 -type f \( -iname "*.fa" -o -iname "*.fasta" -o -iname "*.fas" -o -iname "*.aln" \) -print0)

if [ ${#fasta_files[@]} -eq 0 ]; then
    echo "No FASTA files found in $INPUT_DIR"
    exit 0
fi

for fasta in "${fasta_files[@]}"; do
    basename=$(basename "$fasta")
    name="${basename%.*}"
    subdir="$OUTPUT_DIR/$name"
    mkdir -p "$subdir"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Building tree for $basename ..."
    iqtree \
        -s "$fasta" \
        -m MFP \
        -B 1000 \
        -T "$NPROC" \
        --prefix "$subdir/$name" \
        -redo
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Done: $basename"
    echo
done

echo "All trees built. Results are in $OUTPUT_DIR/"
