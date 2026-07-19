#!/usr/bin/env python3
"""Align every grouped FASTA file with MACSE.

Expected input structure:
    /home/hellstrom/Загрузки/grouped_by_germlines/v/
        group1.fasta
        group2.fasta
        ...

Output structure:
    /home/hellstrom/Документы/MACSE/v/
        group1_aligned.fasta
        group1_aligned_aa.fasta
        group2_aligned.fasta
        ...

The paths are specified directly below; command-line arguments are not needed.
"""

from __future__ import annotations

import subprocess
from pathlib import Path


# Прямые пути: при переносе проекта достаточно изменить эти три строки.
INPUT_DIR = Path("/home/hellstrom/Загрузки/grouped_by_germlines/v")
OUTPUT_DIR = Path("/home/hellstrom/Документы/MACSE/v")
MACSE_BIN = Path("/home/hellstrom/.conda/envs/bio_env/bin/macse")

OUTPUT_SUFFIX = "_aligned"
AA_SUFFIX = "_aa"


def discover_fasta_files(input_dir: Path) -> list[Path]:
    """Return all FASTA files under the input directory."""
    suffixes = {".fasta", ".fa", ".fna", ".ffn", ".faa"}
    files = [
        path
        for path in input_dir.rglob("*")
        if path.is_file() and path.suffix.lower() in suffixes
    ]
    return sorted(files)


def amino_acid_output_path(output_fasta: Path) -> Path:
    """Return the companion amino-acid alignment path for a NT alignment."""
    return output_fasta.with_name(f"{output_fasta.stem}{AA_SUFFIX}{output_fasta.suffix}")


def run_macse(macse_bin: Path, input_fasta: Path, output_fasta: Path) -> None:
    """Run MACSE and write nucleotide and amino-acid alignments."""
    aa_output_fasta = amino_acid_output_path(output_fasta)
    temp_nt_output = output_fasta.with_name(f"{output_fasta.stem}.tmp{output_fasta.suffix}")
    temp_aa_output = aa_output_fasta.with_name(f"{aa_output_fasta.stem}.tmp{aa_output_fasta.suffix}")

    command = [
        str(macse_bin),
        "-prog",
        "alignSequences",
        "-seq",
        str(input_fasta),
        "-out_NT",
        str(temp_nt_output),
        "-out_AA",
        str(temp_aa_output),
    ]

    try:
        subprocess.run(command, check=True)
        temp_nt_output.replace(output_fasta)
        temp_aa_output.replace(aa_output_fasta)
    except Exception:
        # Не оставляем пустые или недописанные файлы после ошибки MACSE.
        output_fasta.unlink(missing_ok=True)
        aa_output_fasta.unlink(missing_ok=True)
        temp_nt_output.unlink(missing_ok=True)
        temp_aa_output.unlink(missing_ok=True)
        raise


def main() -> None:
    input_dir = INPUT_DIR
    output_dir = OUTPUT_DIR

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    if not input_dir.is_dir():
        raise NotADirectoryError(f"Input path is not a directory: {input_dir}")

    if not MACSE_BIN.is_file():
        raise FileNotFoundError(
            f"MACSE not found: {MACSE_BIN}. Install MACSE or change MACSE_BIN "
            "at the beginning of this script."
        )

    fasta_files = discover_fasta_files(input_dir)
    if not fasta_files:
        raise FileNotFoundError(f"No FASTA files found in: {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    manifest_lines = ["input_fasta\taligned_nt_fasta\taligned_aa_fasta"]
    for index, fasta_path in enumerate(fasta_files, start=1):
        rel_path = fasta_path.relative_to(input_dir)
        rel_stem = rel_path.stem
        rel_suffix = rel_path.suffix.lower()
        output_path = output_dir / rel_path.parent / f"{rel_stem}{OUTPUT_SUFFIX}{rel_suffix}"
        output_path.parent.mkdir(parents=True, exist_ok=True)

        aa_output_path = amino_acid_output_path(output_path)

        print(
            f"[{index}/{len(fasta_files)}] Aligning {rel_path} -> "
            f"{output_path.relative_to(output_dir)}"
        )
        run_macse(MACSE_BIN, fasta_path, output_path)
        manifest_lines.append(
            f"{rel_path}\t{output_path.relative_to(output_dir)}\t"
            f"{aa_output_path.relative_to(output_dir)}"
        )

    manifest_path = output_dir / "manifest.tsv"
    manifest_path.write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")

    print()
    print(f"MACSE executable: {MACSE_BIN}")
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Aligned FASTA files: {len(fasta_files)}")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
