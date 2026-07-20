#!/usr/bin/env python3
"""Align FASTA files with MACSE.

Usage:
    python3 MACSEtry.py -i /path/to/input -o /path/to/output
    python3 MACSEtry.py -i /path/to/input -o /path/to/output --macse /opt/macse/macse

MACSE is required (conda install -c bioconda macse or provide --macse).
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


OUTPUT_SUFFIX = "_aligned"


def amino_acid_output_path(output_fasta: Path) -> Path:
    return output_fasta.with_name(f"{output_fasta.stem}_aa{output_fasta.suffix}")


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


def _clean_macse_nt(fasta: Path) -> None:
    text = fasta.read_text(encoding="utf-8")
    if "!" not in text:
        return
    cleaned = text.replace("!", "-")
    fasta.write_text(cleaned, encoding="utf-8")


def run_macse(macse_bin: Path, input_fasta: Path, output_fasta: Path) -> None:
    aa_output = output_fasta.with_name(f"{output_fasta.stem}_aa{output_fasta.suffix}")
    temp_nt = output_fasta.with_name(f"{output_fasta.stem}.tmp{output_fasta.suffix}")
    temp_aa = aa_output.with_name(f"{aa_output.stem}.tmp{aa_output.suffix}")

    cmd = [
        str(macse_bin),
        "-prog", "alignSequences",
        "-seq", str(input_fasta),
        "-out_NT", str(temp_nt),
        "-out_AA", str(temp_aa),
    ]
    try:
        subprocess.run(cmd, check=True)
        temp_nt.replace(output_fasta)
        temp_aa.replace(aa_output)
        _clean_macse_nt(output_fasta)
    except Exception:
        for p in (output_fasta, aa_output, temp_nt, temp_aa):
            p.unlink(missing_ok=True)
        raise


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
    parser = argparse.ArgumentParser(description="Align grouped FASTA files with MACSE.")
    parser.add_argument("-i", "--input", required=True, type=Path,
                        help="Directory with input FASTA files.")
    parser.add_argument("-o", "--output", required=True, type=Path,
                        help="Directory for aligned FASTA files.")
    parser.add_argument("--macse", type=Path, default=None,
                        help="Path to MACSE binary (default: search PATH).")
    parser.add_argument("--threads", type=int, default=1,
                        help="Threads per MACSE instance (default: 1).")
    parser.add_argument("--workers", type=int, default=0,
                        help="Parallel processes (default: CPU count).")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_dir = args.input.expanduser()
    output_dir = args.output.expanduser()
    macse_bin = find_macse(args.macse)
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

    manifest_lines = ["input_fasta\taligned_nt_fasta\taligned_aa_fasta"]
    completed = 0
    total = len(tasks)

    with ThreadPoolExecutor(max_workers=workers) as executor:
        fut_to_task = {
            executor.submit(run_macse, macse_bin, fp, op): (fp, rp, op)
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
            aa_path = op.with_name(f"{op.stem}_aa{op.suffix}")
            print(f"[{completed}/{total}] Aligned {rp}")
            manifest_lines.append(
                f"{rp}\t{op.relative_to(output_dir)}\t{aa_path.relative_to(output_dir)}"
            )

    manifest_path = output_dir / "manifest.tsv"
    manifest_path.write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")

    print()
    print(f"MACSE executable: {macse_bin}")
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Workers: {workers}")
    print(f"Aligned FASTA files: {completed}")
    print(f"Failed: {total - completed}")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
