#!/usr/bin/env python3
"""Batch runner for mutation analysis: format BLAST DBs, then run analyze_mutations.py on each FASTA.

Usage:
    python run_mutations.py -i fasta_from_clades -o mutation_tables -r data/references
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


FASTA_EXTS = {".fa", ".fasta", ".fas", ".aln"}
BLAST_DB_NAMES = ["all_V", "all_D", "all_J"]


def discover_fasta_files(input_dir: Path) -> list[Path]:
    files = [
        p for p in input_dir.iterdir()
        if p.is_file() and p.suffix.lower() in FASTA_EXTS
    ]
    return sorted(files)


def ensure_blast_dbs(ref_dir: Path) -> None:
    for dbname in BLAST_DB_NAMES:
        nsq = ref_dir / f"{dbname}.nsq"
        if not nsq.exists():
            fasta = ref_dir / f"{dbname}.fasta"
            if not fasta.exists():
                print(f"Warning: {fasta} not found, skipping database format", file=sys.stderr)
                continue
            print(f"Formatting BLAST database: {dbname} ...")
            subprocess.run([
                "makeblastdb", "-in", str(fasta),
                "-dbtype", "nucl",
                "-out", str(ref_dir / dbname),
                "-parse_seqids",
            ], check=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Batch mutation analysis: format BLAST DBs and analyze each FASTA."
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Input directory with FASTA files (from clades)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output directory for mutation tables")
    parser.add_argument("-r", "--ref-dir", required=True,
                        help="Reference germline database directory")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input)
    output_dir = Path(args.output)
    ref_dir = Path(args.ref_dir)

    if not input_dir.is_dir():
        print(f"Error: directory {input_dir} not found", file=sys.stderr)
        sys.exit(1)
    if not ref_dir.is_dir():
        print(f"Error: reference directory {ref_dir} not found", file=sys.stderr)
        sys.exit(1)

    if not shutil.which("igblastn"):
        print("Error: igblastn not found. Activate the pipeline environment.", file=sys.stderr)
        sys.exit(1)
    if not shutil.which("makeblastdb"):
        print("Error: makeblastdb not found. Activate the pipeline environment.", file=sys.stderr)
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    ensure_blast_dbs(ref_dir)

    fasta_files = discover_fasta_files(input_dir)
    if not fasta_files:
        print("No FASTA files found in", input_dir)
        sys.exit(0)

    script_dir = Path(__file__).resolve().parent
    analyzer = script_dir / "analyze_mutations.py"

    for fasta in fasta_files:
        name = fasta.stem
        subdir = output_dir / name
        subdir.mkdir(parents=True, exist_ok=True)

        print(f"Analyzing {fasta.name} ...")

        subprocess.run([
            sys.executable, str(analyzer),
            "-i", str(fasta),
            "-o", str(subdir / "mutations.tsv"),
            "--ref-dir", str(ref_dir),
            "--format", "mutations",
        ], check=True)

        subprocess.run([
            sys.executable, str(analyzer),
            "-i", str(fasta),
            "-o", str(subdir / "mutations_summary.tsv"),
            "--ref-dir", str(ref_dir),
            "--format", "summary",
        ], check=True)

        print(f"Done: {fasta.name}")
        print()

    print("All mutation tables saved to", output_dir)


if __name__ == "__main__":
    main()
