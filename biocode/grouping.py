"""Группировка последовательностей в клональные семейства и выбор outgroup.

Группа = общий V(+J)-ген (настраивается ``RunConfig.group_by``). Для каждой группы
строится outgroup = зародышевая (germline) последовательность-предок — это корень
дерева, задающий направление эволюции germline → зрелое антитело.
"""
from __future__ import annotations

import re
from collections import defaultdict

from . import annotate
from .config import RunConfig
from .logging_ import get_logger
from .model import Group, SequenceRecord

log = get_logger("grouping")

_SAFE = re.compile(r"[^A-Za-z0-9._-]+")


def _key_fields(rec: SequenceRecord, group_by: str) -> tuple[str, str]:
    if group_by == "v_gene":
        return rec.v_gene, ""
    if group_by == "v_call":
        return rec.v_call, ""
    if group_by == "v_j_call":
        return rec.v_call, rec.j_call
    return rec.v_gene, rec.j_gene          # v_j_gene (дефолт)


def _sanitize(*parts: str) -> str:
    key = "_".join(p for p in parts if p)
    return _SAFE.sub("-", key) or "UNKNOWN"


def build_outgroup(key: str, records: list[SequenceRecord]) -> SequenceRecord | None:
    """Germline-предок группы как отдельная запись-outgroup (корень дерева)."""
    germ = annotate.germline_consensus([r.germline or "" for r in records])
    if not germ:
        return None
    return SequenceRecord(id=f"GERMLINE_{key}", seq=germ, locus=records[0].locus,
                          productive=True, meta={"synthetic": True, "role": "outgroup"})


def group_records(records: list[SequenceRecord], cfg: RunConfig) -> list[Group]:
    """Разбить записи на группы по ключу; построить outgroup; отсортировать по размеру."""
    buckets: dict[tuple[str, str], list[SequenceRecord]] = defaultdict(list)
    for r in records:
        buckets[_key_fields(r, cfg.group_by)].append(r)

    groups: list[Group] = []
    for (v, j), recs in buckets.items():
        key = _sanitize(v, j)
        groups.append(Group(key=key, v_gene=v, j_gene=j, records=recs,
                            outgroup=build_outgroup(key, recs)))
    groups.sort(key=lambda g: g.size, reverse=True)

    n_small = sum(1 for g in groups if g.size < cfg.min_group_size)
    log.info("групп: %d (из них меньше min_group_size=%d: %d)",
             len(groups), cfg.min_group_size, n_small)
    if cfg.limit_groups > 0:
        groups = groups[:cfg.limit_groups]
    return groups


def analyzable(groups: list[Group], cfg: RunConfig) -> tuple[list[Group], list[Group]]:
    """Разделить на пригодные (size ≥ min) и пропускаемые (с причиной в логе)."""
    ok, skip = [], []
    for g in groups:
        (ok if g.size >= cfg.min_group_size else skip).append(g)
    for g in skip:
        log.debug("skip группа %s: size=%d < %d", g.key, g.size, cfg.min_group_size)
    return ok, skip
