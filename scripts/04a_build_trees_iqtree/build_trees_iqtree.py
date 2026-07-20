#!/usr/bin/env python3
"""Build IQ-TREE trees for aligned FASTA files.

Usage:
    python build_trees_iqtree.py -i aligned_sequences -o trees [--model GTR+F+I+G4]
"""

import argparse
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from shared.utils import count_sequences, discover_fasta_files, detect_nproc

TEMP_IQTREE_EXTS = {".log", ".bionj", ".mldist", ".iqtree", ".ckp.gz",
                    ".model.gz", ".splits.nex", ".vcf", ".contree"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build IQ-TREE trees for aligned FASTA files."
    )
    parser.add_argument("-i", "--input", default="aligned_sequences",
                        help="Input directory with aligned FASTA files")
    parser.add_argument("-o", "--output", default="trees",
                        help="Output directory for tree files")
    parser.add_argument("--model", default="GTR+F+I+G4",
                        help="Substitution model (default: GTR+F+I+G4, use MFP for auto)")
    parser.add_argument("--workers", type=int, default=0,
                        help="Number of parallel IQ-TREE processes (default: nproc // 4)")
    return parser.parse_args()


def run_iqtree(fasta: Path, name: str, subdir: Path, model: str, threads: int) -> None:
    subdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "iqtree",
        "-s", str(fasta),
        "-m", model,
        "-B", "1000",
        "-T", str(threads),
        "--prefix", str(subdir / name),
        "-redo",
    ]
    subprocess.run(cmd, check=True)
    for f in subdir.iterdir():
        if f.suffix in TEMP_IQTREE_EXTS or (f.suffix == ".gz" and not f.name.endswith(".treefile")):
            f.unlink(missing_ok=True)


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input)
    output_dir = Path(args.output)
    nproc = detect_nproc()
    workers = args.workers if args.workers > 0 else max(1, nproc // 4)
    iqtree_threads = max(1, nproc // workers)

    if not input_dir.is_dir():
        print(f"Error: directory {input_dir} not found", file=sys.stderr)
        sys.exit(1)

    if not shutil.which("iqtree"):
        print("Error: iqtree not found. Activate the pipeline environment.", file=sys.stderr)
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    fasta_files = discover_fasta_files(input_dir)
    if not fasta_files:
        print("No FASTA files found in", input_dir)
        sys.exit(0)

    tasks = []
    for fasta in fasta_files:
        name = fasta.stem
        nseqs = count_sequences(fasta)
        if nseqs < 4:
            print(f"[{name}] Skipping: only {nseqs} sequences (need >=4)")
            continue
        subdir = output_dir / name
        tasks.append((fasta, name, subdir))

    total = len(tasks)
    completed = 0

    print(f"Model: {args.model}, workers: {workers}, threads per job: {iqtree_threads}")

    with ThreadPoolExecutor(max_workers=workers) as executor:
        fut_to_task = {
            executor.submit(run_iqtree, fp, nm, sd, args.model, iqtree_threads): (fp, nm, sd)
            for fp, nm, sd in tasks
        }
        for future in as_completed(fut_to_task):
            fp, nm, sd = fut_to_task[future]
            try:
                future.result()
            except Exception as e:
                print(f"[FAILED] {nm}: {e}", file=sys.stderr)
                continue
            completed += 1
            print(f"[{completed}/{total}] Built {nm}")

    print(f"All trees built. Results are in {output_dir}")


if __name__ == "__main__":
    main()
