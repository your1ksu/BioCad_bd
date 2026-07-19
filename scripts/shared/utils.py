from __future__ import annotations

import os
import sys
from pathlib import Path

CODON_TABLE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}


def translate(codon: str) -> str:
    codon = codon.upper().replace('U', 'T')
    return CODON_TABLE.get(codon, 'X')


def count_sequences(fasta_path: Path) -> int:
    count = 0
    with open(fasta_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


FASTA_EXTS = {".fa", ".fasta", ".fas", ".fna", ".ffn", ".faa", ".aln"}


def discover_fasta_files(input_dir: Path, exts: set[str] | None = None) -> list[Path]:
    if exts is None:
        exts = FASTA_EXTS
    return sorted(
        p for p in input_dir.iterdir()
        if p.is_file() and p.suffix.lower() in exts
    )


def detect_nproc() -> int:
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        pass
    n = os.cpu_count()
    return n if n else 1


def parse_gene_name(full_name: str) -> tuple[str, str | None]:
    parts = full_name.split("*")
    gene = parts[0]
    allele = parts[1] if len(parts) > 1 else None
    return gene, allele


def format_group_key(v_name: str, j_name: str, strategy: str) -> str:
    if strategy == "v_only":
        v = v_name.split("*")[0]
        return v
    if strategy == "gene":
        v = v_name.split("*")[0]
        j = j_name.split("*")[0]
        return f"{v}_{j}"
    return f"{v_name}_{j_name}"


def read_germline_fasta(path: str | Path) -> dict[str, str]:
    seqs = {}
    name = None
    chunks = []

    def flush():
        if name is not None:
            seqs[name] = "".join(chunks).upper()

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                parts = line[1:].split("|")
                name = parts[1] if len(parts) > 1 else parts[0]
                chunks = []
            else:
                chunks.append(line)
        flush()
    return seqs


def read_fasta(path: str | Path) -> dict[str, str]:
    records = {}
    key = None
    seq = []

    def flush():
        if key is not None:
            records[key] = "".join(seq)

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                key = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        flush()
    return records


def safe_filename(name: str) -> str:
    import re
    return re.sub(r'[^A-Za-z0-9_\-\.]', '_', name)
