#!/usr/bin/env python3
"""Find SHM hotspots and coldspots in FASTA files and output a JSON report.

Supports gap-aware matching and records precise coordinate boundaries.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path


FASTA_SUFFIXES = {".fasta", ".fa", ".fna", ".ffn"}

IUPAC: dict[str, set[str]] = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "W": {"A", "T"},
    "S": {"C", "G"},
}

HOTSPOT_MOTIFS = ("RGYW", "WRCY", "TAA")
COLDSPOT_MOTIFS = ("SYC", "GGR")


@dataclass(frozen=True)
class FastaRecord:
    header: str
    sequence: str


@dataclass(frozen=True)
class MotifMatch:
    type: str          
    motif_name: str    
    sequence: str      
    start_index: int   
    end_index: int     
    gapped: bool       


def read_fasta(path: Path) -> list[FastaRecord]:
    """Read FASTA records from path."""
    records: list[FastaRecord] = []
    header: str | None = None
    sequence_lines: list[str] = []

    with path.open(encoding="utf-8") as fasta_file:
        for line_number, raw_line in enumerate(fasta_file, start=1):
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append(FastaRecord(header, "".join(sequence_lines)))
                header = line[1:].strip()
                sequence_lines = []
                continue

            if header is None:
                raise ValueError(
                    f"Sequence line before first FASTA header in {path}:{line_number}"
                )
            sequence_lines.append(line)

    if header is not None:
        records.append(FastaRecord(header, "".join(sequence_lines)))

    return records


def motif_matches(window: str, motif: str) -> bool:
    """Return True if sequence window matches an IUPAC motif."""
    if len(window) != len(motif):
        return False

    for nucleotide, motif_symbol in zip(window.upper(), motif):
        allowed = IUPAC.get(motif_symbol, set())
        if nucleotide not in allowed:
            return False

    return True


def find_motifs(sequence: str, motifs: tuple[str, ...], motif_type: str) -> list[MotifMatch]:
    """Find all motif matches, ignoring gap characters (-) for matching logic."""
    matches: list[MotifMatch] = []
    seq_len = len(sequence)

    for motif in motifs:
        motif_length = len(motif)
        
        for start in range(seq_len):
            window_chars: list[str] = []
            current_index = start
            has_gaps = False
            
            # ИСПРАВЛЕНО: Безопасный цикл без риска бесконечного зависания
            while len(window_chars) < motif_length and current_index < seq_len:
                char = sequence[current_index]
                if char == '-':
                    has_gaps = True
                else:
                    window_chars.append(char)
                current_index += 1  # Индекс гарантированно увеличивается всегда
            
            if len(window_chars) == motif_length:
                window_str = "".join(window_chars)
                if motif_matches(window_str, motif):
                    matches.append(
                        MotifMatch(
                            type=motif_type,
                            motif_name=motif,
                            sequence=window_str,
                            start_index=start,
                            end_index=current_index,
                            gapped=has_gaps
                        )
                    )

    return matches


def process_fasta_file(input_path: Path) -> list[dict]:
    """Process a single FASTA file and return list of sequence data dictionaries."""
    records = read_fasta(input_path)
    results = []

    for record in records:
        hotspots = find_motifs(record.sequence, HOTSPOT_MOTIFS, "hotspot")
        coldspots = find_motifs(record.sequence, COLDSPOT_MOTIFS, "coldspot")
        
        all_matches = hotspots + coldspots
        all_matches.sort(key=lambda m: (m.start_index, m.end_index))

        # Берем ID до первого пробела
        seq_id = record.header.split()[0] if record.header.split() else "unknown"

        results.append({
            "sequence_id": seq_id,
            "description": record.header,
            "total_hotspots": len(hotspots),
            "total_coldspots": len(coldspots),
            "matches": [asdict(m) for m in all_matches]
        })

    return results


def discover_fasta_files(input_dir: Path) -> list[Path]:
    """Return all FASTA files under input_dir."""
    return sorted(
        path
        for path in input_dir.rglob("*")
        if path.is_file() and path.suffix.lower() in FASTA_SUFFIXES
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Find SHM hotspots/coldspots and save results to a JSON file."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=Path,
        help="Input FASTA file or directory containing FASTA files.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="Output JSON file path (e.g., report.json).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = args.input.expanduser()
    output_path = args.output.expanduser()

    if not input_path.exists():
        raise FileNotFoundError(f"Input path not found: {input_path}")

    all_results = []

    if input_path.is_file():
        all_results.extend(process_fasta_file(input_path))
    elif input_path.is_dir():
        fasta_files = discover_fasta_files(input_path)
        if not fasta_files:
            raise FileNotFoundError(f"No FASTA files found in: {input_path}")
        for fasta_path in fasta_files:
            all_results.extend(process_fasta_file(fasta_path))
    else:
        raise ValueError(f"Invalid input path type: {input_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as json_file:
        json.dump(all_results, json_file, indent=4, ensure_ascii=False)

    total_seqs = len(all_results)
    total_hot = sum(r["total_hotspots"] for r in all_results)
    total_cold = sum(r["total_coldspots"] for r in all_results)
    
    print(f"Обработано последовательностей: {total_seqs}")
    print(f"Найдено горячих точек (Hotspots): {total_hot}")
    print(f"Найдено холодных точек (Coldspots): {total_cold}")
    print(f"Результаты сохранены в: {output_path}")


if __name__ == "__main__":
    main()
