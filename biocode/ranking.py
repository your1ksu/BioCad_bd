"""Метрики значимости мутаций и ранжирование кандидатов (блок B9).

Ветвевые мутации (:class:`Mutation`) агрегируются по сигнатуре замены
``(position, ref, alt)`` внутри группы и оцениваются композитным скором из
прозрачных признаков (каждый — с весом из ``RunConfig.weights``):

  • recurrence  — независимые повторы одной и той же замены (гомоплазия →
                  сигнал положительного отбора / конвергенции);
  • region_cdr  — регион: CDR3 > CDR1/2 > FR (антиген-связывающие петли важнее);
  • replacement — несинонимичность (замена аминокислоты важнее «молчащей»);
  • support     — надёжность ветви, на которой произошла замена (UFBoot/aLRT);
  • persistence — доля внутренних ветвей (замена наследуется кладой, а не приватна).

Скор — не «чёрный ящик»: возвращаем раскладку вкладов каждого признака.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from .model import Mutation

_CDR = {"CDR1", "CDR2", "CDR3"}
_REGION_WEIGHT = {"CDR3": 1.0, "CDR1": 0.7, "CDR2": 0.7,
                  "FR1": 0.2, "FR2": 0.2, "FR3": 0.2, "FR4": 0.2, "other": 0.0}
_RECURRENCE_CAP = 5


@dataclass
class Candidate:
    """Агрегированная мутация-кандидат со скором и раскладкой."""
    group: str
    position: int
    ref: str
    alt: str
    region: str
    kind: str
    n_branches: int            # recurrence (сколько независимых ветвей)
    n_internal: int            # из них внутренних (persistence)
    max_support: float         # макс. UFBoot по ветвям (0..100)
    branches: list[str]
    score: float = 0.0
    contributions: dict[str, float] = field(default_factory=dict)

    def as_row(self) -> dict:
        r = {k: getattr(self, k) for k in
             ("group", "position", "ref", "alt", "region", "kind",
              "n_branches", "n_internal", "max_support", "score")}
        r["branches"] = ";".join(self.branches)
        r["contributions"] = ";".join(f"{k}={v:.3f}" for k, v in self.contributions.items())
        return r


def aggregate(muts: list[Mutation]) -> list[Candidate]:
    """Свернуть ветвевые мутации в кандидатов по (position, ref, alt)."""
    agg: dict[tuple, dict] = {}
    for m in muts:
        key = (m.position, m.ref, m.alt)
        a = agg.setdefault(key, {
            "group": m.group, "region": m.region, "kind": m.kind,
            "branches": [], "n_internal": 0, "max_support": 0.0})
        a["branches"].append(m.branch)
        if m.branch_internal:
            a["n_internal"] += 1
        a["max_support"] = max(a["max_support"], float(m.support.get("ufboot", 0.0)))
        # регион/тип берём наиболее «сильные» из встреченных
        if _REGION_WEIGHT.get(m.region, 0) > _REGION_WEIGHT.get(a["region"], 0):
            a["region"] = m.region
        if m.kind == "replacement":
            a["kind"] = "replacement"
    out = []
    for (pos, ref, alt), a in agg.items():
        out.append(Candidate(
            group=a["group"], position=pos, ref=ref, alt=alt,
            region=a["region"], kind=a["kind"],
            n_branches=len(a["branches"]), n_internal=a["n_internal"],
            max_support=a["max_support"], branches=a["branches"]))
    return out


def _features(c: Candidate) -> dict[str, float]:
    f_recurrence = min(c.n_branches, _RECURRENCE_CAP) / _RECURRENCE_CAP
    f_region = _REGION_WEIGHT.get(c.region, 0.0)
    f_replacement = 1.0 if c.kind == "replacement" else (0.3 if c.kind == "unknown" else 0.0)
    f_support = c.max_support / 100.0
    f_persistence = (c.n_internal / c.n_branches) if c.n_branches else 0.0
    return {"recurrence": f_recurrence, "region_cdr": f_region,
            "replacement": f_replacement, "support": f_support,
            "persistence": f_persistence}


def score(candidates: list[Candidate], weights: dict[str, float]) -> list[Candidate]:
    """Проставить скор и раскладку вкладов; вернуть отсортированными по убыванию."""
    for c in candidates:
        feats = _features(c)
        contrib = {k: round(weights.get(k, 0.0) * v, 4) for k, v in feats.items()}
        c.contributions = contrib
        c.score = round(sum(contrib.values()), 4)
    candidates.sort(key=lambda c: (c.score, c.n_branches, c.max_support), reverse=True)
    return candidates


def rank_mutations(muts: list[Mutation], weights: dict[str, float]) -> list[Candidate]:
    """Полный конвейер: агрегация → скоринг → сортировка."""
    return score(aggregate(muts), weights)
