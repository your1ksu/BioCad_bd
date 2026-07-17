"""Аннотация FR/CDR и работа с germline.

Каноничный (и точный, без внешних зависимостей) путь для наших AIRR-данных:
столбцы fwr1..fwr4/cdr1..cdr3 конкатенируются побуквенно в ``sequence_alignment``
(проверено), поэтому границы FR/CDR — точные смещения, а IMGT-нумерация делает
FR1–FR3 одинаковыми у всех членов V-группы. Метки «прилеплены» к остаткам и
переносятся через MSA (:func:`project_onto_msa`).

Запасной путь для «голых» FASTA (нет AIRR-разметки) — IgBLAST/ANARCI; обёртки —
в ``structural_platform/predict`` (реализуется в блоке B2, см. PHYLO_PLAN.md).

Standalone-инструмент визуализации FR/CDR (HTML) — ``EDU/region_annotate.py``.
"""
from __future__ import annotations

from collections import Counter

AIRR_REGIONS = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]
LABEL = {"fwr1": "FR1", "cdr1": "CDR1", "fwr2": "FR2", "cdr2": "CDR2",
         "fwr3": "FR3", "cdr3": "CDR3", "fwr4": "FR4"}
REGION_ORDER = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]
CDR_REGIONS = {"CDR1", "CDR2", "CDR3"}
GAP_CHARS = "-."


def _strip_gaps(s: str) -> str:
    return s.replace("-", "").replace(".", "")


def has_airr_regions(row: dict) -> bool:
    return any(row.get(c) for c in AIRR_REGIONS)


def labeled_residues(row: dict) -> tuple[str, list[str]]:
    """Ungapped-последовательность (из region-столбцов) + метка FR/CDR на каждый остаток."""
    seq_chars: list[str] = []
    labels: list[str] = []
    for col in AIRR_REGIONS:
        for ch in (row.get(col) or ""):
            if ch in GAP_CHARS:
                continue
            seq_chars.append(ch.upper())
            labels.append(LABEL[col])
    return "".join(seq_chars), labels


def regions_from_row(row: dict) -> dict[str, tuple[int, int]] | None:
    """Границы FR/CDR как {'FR1': (start,end), ...} в координатах UNGAPPED-остатков.

    Валидно, когда конкатенация region-столбцов == sequence_alignment (наши данные).
    Возвращает None, если разметки нет.
    """
    if not has_airr_regions(row):
        return None
    spans: dict[str, tuple[int, int]] = {}
    pos = 0
    for col in AIRR_REGIONS:
        piece = _strip_gaps(row.get(col) or "")
        if piece:
            spans[LABEL[col]] = (pos, pos + len(piece))
            pos += len(piece)
    return spans or None


def regions_valid(row: dict) -> bool:
    """True, если region-столбцы побуквенно = sequence_alignment (разметке можно доверять)."""
    concat = "".join((row.get(c) or "") for c in AIRR_REGIONS)
    return bool(concat) and concat == (row.get("sequence_alignment") or "")


def project_onto_msa(labeled: dict[str, list[str]], msa: dict[str, str]) -> list[dict]:
    """Спроецировать FR/CDR-метки на колонки готового MSA.

    labeled: id → метки на каждый негэповый остаток (из :func:`labeled_residues`).
    msa:     id → строка выравнивания (с '-'). Метки едут с остатками → проекция точна.
    Возвращает по колонке: {column, region (доминирующий), purity, coverage}.
    """
    width = max((len(s) for s in msa.values()), default=0)
    cols = [Counter() for _ in range(width)]
    for sid, aligned in msa.items():
        labs = labeled.get(sid)
        if not labs:
            continue
        k = 0
        for c, ch in enumerate(aligned):
            if ch in GAP_CHARS:
                continue
            if k < len(labs):
                cols[c][labs[k]] += 1
            k += 1
    out = []
    for c, cnt in enumerate(cols):
        total = sum(cnt.values())
        if total == 0:
            out.append({"column": c, "region": "other", "purity": 0.0, "coverage": 0})
        else:
            region, n = cnt.most_common(1)[0]
            out.append({"column": c, "region": region,
                        "purity": round(n / total, 3), "coverage": total})
    return out


def region_at_column(track: list[dict], column: int) -> str:
    """Регион для колонки MSA по треку project_onto_msa (для разметки мутаций)."""
    if 0 <= column < len(track):
        return track[column]["region"]
    return "other"


def germline_consensus(germlines: list[str]) -> str | None:
    """Представительный germline группы = самый частый непустой germline-стринг.

    В пределах одной V+J-группы зародышевые последовательности почти идентичны,
    поэтому мода — корректный предок (корень дерева). Fallback — самый длинный.
    """
    cands = [g for g in germlines if g]
    if not cands:
        return None
    counts = Counter(cands)
    top, n = counts.most_common(1)[0]
    if n > 1:
        return top
    return max(cands, key=len)
