"""Сравнение ML- и байесовского деревьев (блок B7).

Метрики согласия двух деревьев одной группы:
  • Robinson–Foulds — число несовпадающих биразбиений (топологическая дистанция);
  • clade agreement — общие биразбиения + опоры на них (UFBoot/aLRT vs posterior);
  • posterior_by_signature — апостериорные по подписи клады (для «уверенных в обеих
    моделях» в clades.py).

Реализация на биразбиениях через ``Bio.Phylo`` (без dendropy/ete3).
"""
from __future__ import annotations

from io import StringIO

from ..model import TreeResult


def _leaves(tree) -> set[str]:
    return {t.name for t in tree.get_terminals()}


def _splits(newick: str, exclude: set[str] | None = None) -> tuple[set[frozenset], set[str]]:
    """Множество нетривиальных биразбиений дерева.

    Каждое внутреннее ребро делит листья на две стороны; берём сторону БЕЗ якоря
    (детерминированно), чтобы одно и то же разбиение имело одинаковое представление
    в обоих деревьях. Якорь — первый по сортировке общий лист.
    """
    from Bio import Phylo

    tree = Phylo.read(StringIO(newick), "newick")
    exclude = exclude or set()
    all_leaves = _leaves(tree) - exclude
    if not all_leaves:
        return set(), set()
    anchor = sorted(all_leaves)[0]
    splits: set[frozenset] = set()
    n = len(all_leaves)
    for clade in tree.get_nonterminals():
        side = frozenset(t.name for t in clade.get_terminals()) - exclude
        if len(side) < 2 or len(side) > n - 2:
            continue
        canon = side if anchor not in side else (all_leaves - side)
        splits.add(frozenset(canon))
    return splits, all_leaves


def robinson_foulds(ml: TreeResult, bayes: TreeResult) -> dict:
    """RF-дистанция между ML- и байесовским деревьями (без outgroup)."""
    exclude = {ml.outgroup, bayes.outgroup} - {None}
    s1, l1 = _splits(ml.newick, exclude)
    s2, l2 = _splits(bayes.newick, exclude)
    if l1 != l2:
        # сравниваем только по общему набору листьев
        common = l1 & l2
        s1 = {frozenset(s & common) for s in s1 if len(s & common) >= 2}
        s2 = {frozenset(s & common) for s in s2 if len(s & common) >= 2}
        n = len(common)
    else:
        n = len(l1)
    rf = len(s1 ^ s2)
    max_rf = 2 * (n - 3) if n >= 3 else 0
    return {
        "rf": rf,
        "max_rf": max_rf,
        "normalized_rf": round(rf / max_rf, 4) if max_rf else 0.0,
        "shared_splits": len(s1 & s2),
        "ml_only": len(s1 - s2),
        "bayes_only": len(s2 - s1),
        "n_taxa": n,
    }


def posterior_by_signature(bayes: TreeResult) -> dict[str, float]:
    """{clade_sig → posterior} из supports байесовского дерева (для clades.py)."""
    return {sig: s["posterior"] for sig, s in bayes.supports.items()
            if "posterior" in s}


def clade_agreement(ml: TreeResult, bayes: TreeResult) -> list[dict]:
    """Общие клады: их поддержка в ML (ufboot/alrt) и апостериорная в Bayes."""
    post = posterior_by_signature(bayes)
    rows = []
    for sig, ml_sup in ml.supports.items():
        if sig in post:
            rows.append({
                "clade_size": sig.count("|") + 1,
                "ufboot": ml_sup.get("ufboot"),
                "alrt": ml_sup.get("alrt"),
                "posterior": post[sig],
            })
    rows.sort(key=lambda r: (r["posterior"], r.get("ufboot") or 0), reverse=True)
    return rows


def compare_trees(ml: TreeResult, bayes: TreeResult) -> dict:
    """Полная сводка согласия ML↔Bayes для отчёта группы."""
    agree = clade_agreement(ml, bayes)
    return {
        "robinson_foulds": robinson_foulds(ml, bayes),
        "shared_clades": len(agree),
        "shared_clade_support": agree,
    }
