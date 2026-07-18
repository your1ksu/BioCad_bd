"""Ввод/вывод: AIRR-TSV и FASTA → SequenceRecord[]; атомарная запись артефактов."""
from __future__ import annotations

import csv
import json
import os
import tempfile
from pathlib import Path

from . import annotate
from .config import RunConfig
from .errors import InputError
from .logging_ import get_logger
from .model import SequenceRecord

log = get_logger("io")

_TRUE = {"t", "true", "1", "yes"}


def _bool(v) -> bool:
    return str(v).strip().lower() in _TRUE


def _isotype(c_gene: str | None) -> str | None:
    if not c_gene:
        return None
    for iso in ("IGHM", "IGHD", "IGHG", "IGHA", "IGHE", "IGKC", "IGLC"):
        if c_gene.upper().startswith(iso):
            return iso
    return c_gene


def _working_seq(row: dict, working: str) -> str:
    raw = row.get(working) or ""
    return annotate._strip_gaps(raw).upper()


def read_airr_tsv(path: str | Path, cfg: RunConfig) -> list[SequenceRecord]:
    """Загрузить AIRR-TSV в SequenceRecord[]. Фильтры locus/productive — по конфигу."""
    path = Path(path)
    if not path.is_file():
        raise InputError(f"AIRR-TSV не найден: {path}")
    records: list[SequenceRecord] = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None or "sequence_id" not in reader.fieldnames:
            raise InputError(f"{path}: не похоже на AIRR-TSV (нет столбца sequence_id)")
        for row in reader:
            if cfg.locus and (row.get("locus") or "") != cfg.locus:
                continue
            if cfg.productive_only and not _bool(row.get("productive")):
                continue
            seq = _working_seq(row, cfg.working_seq)
            if not seq:
                continue
            # germline для outgroup: d_mask (D-регион замаскирован N) — стандарт для
            # укоренения BCR-линий; fallback на germline_alignment. "NA"/только-N → нет.
            germ_raw = row.get("germline_alignment_d_mask") or row.get("germline_alignment") or ""
            germ = annotate._strip_gaps(germ_raw).upper()
            if germ in ("", "NA", "N/A") or set(germ) <= {"N"}:
                germ = None
            records.append(SequenceRecord(
                id=row["sequence_id"],
                seq=seq,
                v_gene=row.get("v_gene") or "", j_gene=row.get("j_gene") or "",
                d_gene=(row.get("d_gene") or None),
                v_call=row.get("v_call") or "", j_call=row.get("j_call") or "",
                d_call=(row.get("d_call") or None),
                locus=row.get("locus") or "",
                productive=_bool(row.get("productive")),
                isotype=_isotype(row.get("c_gene")),
                regions=annotate.regions_from_row(row),
                germline=germ,
                meta={
                    "cell_id": row.get("cell_id"),
                    "duplicate_count": row.get("duplicate_count"),
                    "day": row.get("day"),
                    "sample_id": row.get("sample_id"),
                    "subject_id": row.get("subject_id"),
                    "regions_valid": annotate.regions_valid(row),
                    "labels": annotate.labeled_residues(row)[1] if annotate.has_airr_regions(row) else None,
                },
            ))
    log.info("загружено %d записей из %s (locus=%s, productive_only=%s)",
             len(records), path.name, cfg.locus or "любой", cfg.productive_only)
    return records


def read_fasta(path: str | Path) -> dict[str, str]:
    """FASTA → {id: seq} (ungapped, upper). id = первое слово заголовка."""
    path = Path(path)
    if not path.is_file():
        raise InputError(f"FASTA не найден: {path}")
    seqs: dict[str, str] = {}
    cur = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            cur = line[1:].strip().split()[0]
            seqs[cur] = ""
        elif cur is not None:
            seqs[cur] += line.strip()
    return {k: v.upper() for k, v in seqs.items()}


def read_fasta_ids(path: str | Path) -> list[str]:
    return list(read_fasta(path).keys())


def records_from_fasta_join(fasta: str | Path, tsv: str | Path,
                            cfg: RunConfig) -> list[SequenceRecord]:
    """Записи по FASTA (напр. группа anotherpipeline), метаданные — join к AIRR-TSV по id."""
    ids = set(read_fasta_ids(fasta))
    all_recs = read_airr_tsv(tsv, RunConfig(**{**cfg.to_dict(),
                                              "locus": "", "productive_only": False}))
    by_id = {r.id: r for r in all_recs}
    missing = ids - set(by_id)
    if missing:
        log.warning("%d/%d id из FASTA не найдены в AIRR-TSV (будут пропущены)",
                    len(missing), len(ids))
    return [by_id[i] for i in ids if i in by_id]


# ─────────────────────────── запись артефактов ────────────────────────────

def _atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=str(path.parent), suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            f.write(text)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def write_fasta(records: list[SequenceRecord], path: str | Path, *, wrap: int = 0) -> Path:
    path = Path(path)
    lines = []
    for r in records:
        lines.append(f">{r.id}")
        s = r.seq
        if wrap and wrap > 0:
            lines += [s[i:i + wrap] for i in range(0, len(s), wrap)]
        else:
            lines.append(s)
    _atomic_write(path, "\n".join(lines) + "\n")
    return path


def write_tsv(rows: list[dict], path: str | Path, columns: list[str] | None = None) -> Path:
    path = Path(path)
    if not rows:
        _atomic_write(path, "\t".join(columns or []) + "\n")
        return path
    cols = columns or list(rows[0].keys())
    out = ["\t".join(cols)]
    for r in rows:
        out.append("\t".join(str(r.get(c, "")) for c in cols))
    _atomic_write(path, "\n".join(out) + "\n")
    return path


def write_json(obj, path: str | Path) -> Path:
    path = Path(path)
    _atomic_write(path, json.dumps(obj, ensure_ascii=False, indent=2, default=str))
    return path
