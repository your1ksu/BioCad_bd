#!/usr/bin/env python3
"""Align every grouped FASTA file with MAFFT.

Example:
    python3 MSA_final.py \
        --input /path/to/grouped_by_germlines \
        --output /path/to/aligned_sequences

If --mafft is not specified, the script searches for "mafft" in PATH.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path


OUTPUT_SUFFIX = "_aligned"


def discover_fasta_files(input_dir: Path) -> list[Path]:
    """Return all FASTA files under the input directory."""
    suffixes = {".fasta", ".fa", ".fna", ".ffn", ".faa"}
    files = [
        path
        for path in input_dir.rglob("*")
        if path.is_file() and path.suffix.lower() in suffixes
    ]
    return sorted(files)


def run_mafft(mafft_bin: Path, input_fasta: Path, output_fasta: Path) -> None:
    """Run MAFFT and write the alignment to the output file."""
    # --preservecase не позволяет MAFFT превращать входные A/C/G/T в строчные.
    command = [
        str(mafft_bin),
        "--auto",
        "--quiet",
        "--preservecase",
        str(input_fasta),
    ]
    try:
        with output_fasta.open("w", encoding="utf-8") as output_handle:
            subprocess.run(command, stdout=output_handle, check=True)
    except Exception:
        # Не оставляем пустой или недописанный файл после ошибки MAFFT.
        output_fasta.unlink(missing_ok=True)
        raise


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Find FASTA files in an input directory and align each file with MAFFT."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=Path,
        help="Directory with input FASTA files.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="Directory where aligned FASTA files and manifest.tsv will be written.",
    )
    parser.add_argument(
        "-m",
        "--mafft",
        type=Path,
        default=None,
        help=(
            "Path to MAFFT executable. If omitted, the script searches for mafft "
            "in PATH."
        ),
    )
    return parser.parse_args()


def find_mafft(mafft_arg: Path | None) -> Path:
    """Return a MAFFT executable from --mafft or from PATH."""
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


def main() -> None:
    args = parse_args()
    input_dir = args.input.expanduser()
    output_dir = args.output.expanduser()
    mafft_bin = find_mafft(args.mafft)

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    if not input_dir.is_dir():
        raise NotADirectoryError(f"Input path is not a directory: {input_dir}")

    fasta_files = discover_fasta_files(input_dir)
    if not fasta_files:
        raise FileNotFoundError(f"No FASTA files found in: {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    manifest_lines = ["input_fasta\taligned_fasta"]
    for index, fasta_path in enumerate(fasta_files, start=1):
        rel_path = fasta_path.relative_to(input_dir)
        rel_stem = rel_path.stem
        rel_suffix = rel_path.suffix.lower()
        output_path = output_dir / rel_path.parent / f"{rel_stem}{OUTPUT_SUFFIX}{rel_suffix}"
        output_path.parent.mkdir(parents=True, exist_ok=True)

        print(f"[{index}/{len(fasta_files)}] Aligning {rel_path} -> {output_path.relative_to(output_dir)}")
        run_mafft(mafft_bin, fasta_path, output_path)
        manifest_lines.append(f"{rel_path}\t{output_path.relative_to(output_dir)}")

    manifest_path = output_dir / "manifest.tsv"
    manifest_path.write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")

    print()
    print(f"MAFFT executable: {mafft_bin}")
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Aligned FASTA files: {len(fasta_files)}")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
