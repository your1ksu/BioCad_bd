"""Уверенные клады — «микроскрипт поиска клад, в которых мы больше всего уверены»
(задача Miro, блок B10).

Из ML-дерева берём внутренние ветви с высокой поддержкой (UFBoot ≥ порог И
aLRT ≥ порог) — это клады, которым можно доверять. Для каждой: состав листьев,
размер, определяющие мутации (замены на входящей в кладу ветви), изотипы/дни
членов. Хук для согласия с Байесом (posterior) — при наличии MrBayes-дерева
(блок B7): клада «уверенная в обеих моделях», если ещё и posterior ≥ 0.95.
"""
from __future__ import annotations

from collections import Counter
from io import StringIO

from .model import Mutation, SequenceRecord, TreeResult


def _node_name(clade) -> str | None:
    if clade.is_terminal():
        return clade.name
    if clade.name:
        return clade.name.split("/")[0]
    return None


def _support(clade) -> dict[str, float]:
    if clade.is_terminal() or not clade.name:
        return {}
    nums = []
    for tok in clade.name.split("/"):
        try:
            nums.append(float(tok))
        except ValueError:
            pass
    if len(nums) >= 2:
        return {"alrt": nums[-2], "ufboot": nums[-1]}
    if len(nums) == 1:
        return {"ufboot": nums[0]}
    return {}


def confident_clades(tree: TreeResult, muts: list[Mutation],
                     records: dict[str, SequenceRecord] | None = None,
                     *, ufboot_min: float = 95.0, alrt_min: float = 80.0,
                     posterior: dict[str, float] | None = None,
                     posterior_min: float = 0.95) -> list[dict]:
    """Список уверенных клад ML-дерева (+ метаданные и определяющие мутации)."""
    from Bio import Phylo

    if not tree.newick:
        return []
    phylo = Phylo.read(StringIO(tree.newick), "newick")
    records = records or {}
    outgroup = tree.outgroup

    by_branch_child: dict[str, list[Mutation]] = {}
    for m in muts:
        child = m.branch.split("→")[-1]
        by_branch_child.setdefault(child, []).append(m)

    clades: list[dict] = []
    for clade in phylo.get_nonterminals():
        node = _node_name(clade)
        sup = _support(clade)
        uf, al = sup.get("ufboot", 0.0), sup.get("alrt", 0.0)
        if not (uf >= ufboot_min and al >= alrt_min):
            continue
        leaves = [t.name for t in clade.get_terminals() if t.name != outgroup]
        if len(leaves) < 2:
            continue
        defin = by_branch_child.get(node or "", [])
        isotypes = Counter(records[i].isotype for i in leaves
                           if i in records and records[i].isotype)
        days = sorted({records[i].meta.get("day") for i in leaves
                       if i in records and records[i].meta.get("day")})
        entry = {
            "clade": node, "size": len(leaves), "leaves": leaves,
            "ufboot": uf, "alrt": al,
            "defining_mutations": len(defin),
            "defining_cdr": sum(1 for m in defin if m.region in ("CDR1", "CDR2", "CDR3")),
            "isotypes": dict(isotypes), "days": days,
            "confident_both_models": False,
        }
        if posterior is not None:
            sig = "|".join(sorted(leaves))
            pp = posterior.get(sig)
            entry["posterior"] = pp
            entry["confident_both_models"] = pp is not None and pp >= posterior_min
        clades.append(entry)
    clades.sort(key=lambda c: (c["ufboot"], c["size"]), reverse=True)
    return clades
