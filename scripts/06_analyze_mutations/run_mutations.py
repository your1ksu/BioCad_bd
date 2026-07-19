#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


FASTA_EXTS = {".fa", ".fasta", ".fas", ".aln"}
BLAST_DB_NAMES = ["all_V", "all_D", "all_J"]


def detect_nproc() -> int:
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        pass
    n = os.cpu_count()
    return n if n else 1


def discover_fasta_files(input_dir: Path) -> list[Path]:
    return sorted(
        p for p in input_dir.iterdir()
        if p.is_file() and p.suffix.lower() in FASTA_EXTS
    )


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


def analyze_one(fasta: Path, analyzer: Path, output_dir: Path, ref_dir: Path) -> str:
    name = fasta.stem
    subdir = output_dir / name
    subdir.mkdir(parents=True, exist_ok=True)

    for fmt, outname in [("mutations", "mutations.tsv"), ("summary", "mutations_summary.tsv")]:
        subprocess.run([
            sys.executable, str(analyzer),
            "-i", str(fasta),
            "-o", str(subdir / outname),
            "--ref-dir", str(ref_dir),
            "--format", fmt,
        ], check=True)

    return name


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
    parser.add_argument("--workers", type=int, default=0,
                        help="Сколько файлов анализировать параллельно (по умолчанию: min(nproc, 4))")
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

    nproc = detect_nproc()
    workers = args.workers if args.workers > 0 else nproc
    workers = min(workers, len(fasta_files))

    script_dir = Path(__file__).resolve().parent
    analyzer = script_dir / "analyze_mutations.py"

    print(f"Analyzing {len(fasta_files)} files with {workers} workers ...")

    completed = 0
    if workers <= 1:
        for fasta in fasta_files:
            name = analyze_one(fasta, analyzer, output_dir, ref_dir)
            completed += 1
            print(f"  [{completed}/{len(fasta_files)}] {name}")
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            fut_to_fasta = {
                pool.submit(analyze_one, fasta, analyzer, output_dir, ref_dir): fasta
                for fasta in fasta_files
            }
            for future in as_completed(fut_to_fasta):
                fasta = fut_to_fasta[future]
                try:
                    name = future.result()
                except Exception as e:
                    print(f"  [FAILED] {fasta.name}: {e}", file=sys.stderr)
                    continue
                completed += 1
                print(f"  [{completed}/{len(fasta_files)}] {name}")

    print(f"\nAll mutation tables saved to {output_dir}")


if __name__ == "__main__":
    main()
