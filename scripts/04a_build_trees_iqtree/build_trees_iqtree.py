#!/usr/bin/env python3
"""Build IQ-TREE trees for aligned FASTA files.

Usage:
    python build_trees_iqtree.py -i aligned_sequences -o trees
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


FASTA_EXTS = {".fa", ".fasta", ".fas", ".aln"}


def discover_fasta_files(input_dir: Path) -> list[Path]:
    files = [
        p for p in input_dir.iterdir()
        if p.is_file() and p.suffix.lower() in FASTA_EXTS
    ]
    return sorted(files)


def count_sequences(fasta_path: Path) -> int:
    count = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def detect_nproc() -> int:
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        pass
    n = os.cpu_count()
    return n if n else 1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build IQ-TREE trees for aligned FASTA files."
    )
    parser.add_argument("-i", "--input", default="aligned_sequences",
                        help="Input directory with aligned FASTA files")
    parser.add_argument("-o", "--output", default="trees",
                        help="Output directory for tree files")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input)
    output_dir = Path(args.output)

    if not input_dir.is_dir():
        print(f"Error: directory {input_dir} not found", file=sys.stderr)
        sys.exit(1)

    if not shutil.which("iqtree"):
        print("Error: iqtree not found. Install it or activate the pipeline environment.", file=sys.stderr)
        sys.exit(1)

    nproc = detect_nproc()
    output_dir.mkdir(parents=True, exist_ok=True)

    fasta_files = discover_fasta_files(input_dir)
    if not fasta_files:
        print("No FASTA files found in", input_dir)
        sys.exit(0)

    for fasta in fasta_files:
        name = fasta.stem
        subdir = output_dir / name
        subdir.mkdir(parents=True, exist_ok=True)

        nseqs = count_sequences(fasta)
        if nseqs < 4:
            print(f"[{name}] Skipping: only {nseqs} sequences (need >=4)")
            continue

        print(f"[{name}] Building tree ({nseqs} seqs) ...")
        cmd = [
            "iqtree",
            "-s", str(fasta),
            "-m", "MFP",
            "-B", "1000",
            "-T", str(nproc),
            "--prefix", str(subdir / name),
            "-redo",
        ]
        subprocess.run(cmd, check=True)
        print(f"[{name}] Done")

    print("All trees built. Results are in", output_dir)


if __name__ == "__main__":
    main()
