#!/usr/bin/env bash

set -euo pipefail

INPUT_DIR="${1:-fasta_from_clades}"
OUTPUT_DIR="${2:-mutation_tables}"
REF_DIR="${3:-}"
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

if ! command -v igblastn &>/dev/null; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Installing igblast, blast, biopython in '${CONDA_ENV}' ..."
    conda install -c conda-forge -c bioconda igblast blast biopython -y
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: directory $INPUT_DIR not found"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

if [ -z "$REF_DIR" ]; then
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REF_DIR="$script_dir/references"
fi

# Format reference databases if needed
for dbname in all_V all_D all_J; do
    nsq="$REF_DIR/${dbname}.nsq"
    if [ ! -f "$nsq" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Formatting BLAST database: $dbname ..."
        makeblastdb -in "$REF_DIR/${dbname}.fasta" -dbtype nucl -out "$REF_DIR/$dbname" -parse_seqids
    fi
done

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

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Analyzing $basename ..."

    python3 "$(dirname "$0")/analyze_mutations.py" \
        -i "$fasta" \
        -o "$subdir/mutations.tsv" \
        --ref-dir "$REF_DIR" \
        --format mutations

    python3 "$(dirname "$0")/analyze_mutations.py" \
        -i "$fasta" \
        -o "$subdir/mutations_summary.tsv" \
        --ref-dir "$REF_DIR" \
        --format summary

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Done: $basename"
    echo
done

echo "All mutation tables saved to $OUTPUT_DIR/"
