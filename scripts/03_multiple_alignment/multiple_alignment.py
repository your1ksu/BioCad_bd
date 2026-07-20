#!/usr/bin/env python3
"""Align every grouped FASTA file with MAFFT or MACSE.

Examples:
    python3 multiple_alignment.py -i /path/to/input -o /path/to/output
    python3 multiple_alignment.py -i /path/to/input -o /path/to/output --aligner macse
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# Allow import of MACSE runner when aligner=macse
import sys
SCRIPTS_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SCRIPTS_DIR))

from MACSE.MACSEtry import run_macse as _run_macse


OUTPUT_SUFFIX = "_aligned"


def detect_nproc() -> int:
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        pass
    n = os.cpu_count()
    return n if n else 1


def discover_fasta_files(input_dir: Path) -> list[Path]:
    suffixes = {".fasta", ".fa", ".fna", ".ffn", ".faa"}
    files = [
        path
        for path in input_dir.rglob("*")
        if path.is_file() and path.suffix.lower() in suffixes
    ]
    return sorted(files)


def run_mafft(mafft_bin: Path, input_fasta: Path, output_fasta: Path, threads: int) -> None:
    command = [
        str(mafft_bin),
        "--auto",
        "--quiet",
        "--preservecase",
        "--thread", str(threads),
        str(input_fasta),
    ]
    try:
        with output_fasta.open("w", encoding="utf-8") as output_handle:
            subprocess.run(command, stdout=output_handle, check=True)
    except Exception:
        output_fasta.unlink(missing_ok=True)
        raise


def run_macse(macse_bin: Path, input_fasta: Path, output_fasta: Path) -> None:
    _run_macse(macse_bin, input_fasta, output_fasta)


def find_mafft(mafft_arg: Path | None) -> Path:
    if mafft_arg is not None:
        mafft_bin = mafft_arg.expanduser()
        if not mafft_bin.is_file():
            raise FileNotFoundError(f"MAFFT not found: {mafft_bin}")
        return mafft_bin
    mafft_from_path = shutil.which("mafft")
    if mafft_from_path is None:
        raise FileNotFoundError(
            "MAFFT not found in PATH. Install MAFFT, activate the environment "
            "that contains it, or pass the executable with --mafft."
        )
    return Path(mafft_from_path)


def find_macse(macse_arg: Path | None) -> Path:
    if macse_arg is not None:
        macse_bin = macse_arg.expanduser()
        if not macse_bin.is_file():
            raise FileNotFoundError(f"MACSE not found: {macse_bin}")
        return macse_bin
    macse_from_path = shutil.which("macse")
    if macse_from_path is None:
        raise FileNotFoundError(
            "MACSE not found in PATH. Install MACSE (conda install -c bioconda macse) "
            "or pass the executable with --macse."
        )
    return Path(macse_from_path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Find FASTA files in an input directory and align each file."
    )
    parser.add_argument("-i", "--input", required=True, type=Path,
                        help="Directory with input FASTA files.")
    parser.add_argument("-o", "--output", required=True, type=Path,
                        help="Directory where aligned FASTA files and manifest.tsv will be written.")
    parser.add_argument("--aligner", choices=["mafft", "macse"], default="mafft",
                        help="Aligner to use (default: mafft).")
    parser.add_argument("-m", "--mafft", type=Path, default=None,
                        help="Path to MAFFT executable (default: search PATH).")
    parser.add_argument("--macse", type=Path, default=None,
                        help="Path to MACSE executable (default: search PATH).")
    parser.add_argument("--threads", type=int, default=1,
                        help="Threads per aligner process (default: 1).")
    parser.add_argument("--workers", type=int, default=0,
                        help="Number of parallel processes (default: CPU count).")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_dir = args.input.expanduser()
    output_dir = args.output.expanduser()
    aligner = args.aligner
    nproc = detect_nproc()
    workers = args.workers if args.workers > 0 else nproc

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    if not input_dir.is_dir():
        raise NotADirectoryError(f"Input path is not a directory: {input_dir}")

    fasta_files = discover_fasta_files(input_dir)
    if not fasta_files:
        raise FileNotFoundError(f"No FASTA files found in: {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    tasks = []
    for fasta_path in fasta_files:
        rel_path = fasta_path.relative_to(input_dir)
        rel_stem = rel_path.stem
        rel_suffix = rel_path.suffix.lower()
        output_path = output_dir / rel_path.parent / f"{rel_stem}{OUTPUT_SUFFIX}{rel_suffix}"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        tasks.append((fasta_path, rel_path, output_path))

    r_mafft = None
    r_macse = None
    if aligner == "mafft":
        r_mafft = find_mafft(args.mafft)
    else:
        r_macse = find_macse(args.macse)

    manifest_lines = ["input_fasta\taligned_fasta"]
    completed = 0
    total = len(tasks)

    with ThreadPoolExecutor(max_workers=workers) as executor:
        if aligner == "mafft":
            fut_to_task = {
                executor.submit(run_mafft, r_mafft, fp, op, args.threads): (fp, rp, op)
                for fp, rp, op in tasks
            }
        else:
            fut_to_task = {
                executor.submit(run_macse, r_macse, fp, op): (fp, rp, op)
                for fp, rp, op in tasks
            }
        for future in as_completed(fut_to_task):
            fp, rp, op = fut_to_task[future]
            try:
                future.result()
            except Exception as e:
                print(f"[FAILED] {rp}: {e}")
                continue
            completed += 1
            print(f"[{completed}/{total}] Aligned {rp}")
            manifest_lines.append(f"{rp}\t{op.relative_to(output_dir)}")

    manifest_path = output_dir / "manifest.tsv"
    manifest_path.write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")

    aligner_name = "MAFFT" if aligner == "mafft" else "MACSE"
    print()
    print(f"{aligner_name} executable: {r_mafft or r_macse}")
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Workers: {workers}")
    print(f"Aligned FASTA files: {completed}")
    print(f"Failed: {total - completed}")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
