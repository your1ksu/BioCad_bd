#!/usr/bin/env python3
"""ШАГ «Никиты» (уверенные клады) — nexus → report.json.

Обёртка над продакшн-движком ``biocode`` (см. ../biocode/, вендорено без
изменений из BIOCAD.bigchallenges@main), адаптированная под файловый контракт
остальных участников конвейера (см. ../BioCad_repo.md, строка «уверенные
клады»: критерий UFBoot≥95 и aLRT≥80, у Байеса — posterior≥0.95).

Два независимых источника NEXUS-деревьев поддержаны (можно использовать один
или оба сразу — тогда клады группы объединяются в один report.json):

  --mrbayes-dir DIR (по умолчанию 'mrbayes')
      <группа>.nex.con.tre из шага ../mrbayes/run_mrbayes.py — консенсус
      MrBayes с апостериорной поддержкой (posterior) на внутренних узлах.
      Критерий: posterior ≥ --posterior-min (по умолчанию 0.95).
      Требуется <группа>.names.tsv рядом (тот же шаг его пишет), чтобы вернуть
      исходные id таксонов вместо MrBayes-совместимых T0001...

  --iqtree-dir DIR (напр. 'trees', выход anotherpipeline/build_trees/build_trees.sh
      у Дениса: trees/<группа>/<группа>.treefile)
      ML-дерево IQ-TREE с метками узлов вида 'SH-aLRT/UFboot' (напр. '98.5/100').
      Критерий: UFBoot ≥ --ufboot-min (95) И aLRT ≥ --alrt-min (80) — здесь
      напрямую вызывается ``biocode.clades.confident_clades`` без изменений.

Выход: groups/report.json — {<группа>: {"mrbayes": {...}, "iqtree": {...}}}
(секции присутствуют только для тех источников, что были переданы и найдены).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from biocode import clades as clades_mod
from biocode.model import TreeResult
from biocode.trees import mrbayes as mrbayes_mod

CON_TRE_SUFFIX = ".nex.con.tre"


def _load_names(names_tsv: Path) -> dict[str, str]:
    """<группа>.names.tsv (safe_id\\toriginal_id) → {safe_id: original_id}."""
    if not names_tsv.is_file():
        return {}
    rev: dict[str, str] = {}
    for line in names_tsv.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        safe, orig = line.split("\t", 1)
        rev[safe] = orig
    return rev


def clades_from_mrbayes(con_tre: Path, names_tsv: Path, posterior_min: float) -> list[dict]:
    """<группа>.nex.con.tre (posterior) → список уверенных клад."""
    rev = _load_names(names_tsv)
    _, supports = mrbayes_mod.parse_con_tre(con_tre, rev)
    entries = []
    for sig, sup in supports.items():
        pp = sup.get("posterior")
        if pp is None or pp < posterior_min:
            continue
        leaves = sig.split("|")
        entries.append({
            "clade": sig, "size": len(leaves), "leaves": leaves,
            "ufboot": None, "alrt": None, "posterior": pp,
            "defining_mutations": None, "defining_cdr": None,
            "isotypes": {}, "days": [], "confident_both_models": False,
        })
    entries.sort(key=lambda c: (c["posterior"], c["size"]), reverse=True)
    return entries


def clades_from_iqtree(treefile: Path, ufboot_min: float, alrt_min: float) -> list[dict]:
    """<группа>.treefile (newick, alrt/ufboot в метках узлов) → уверенные клады.

    Прямой вызов biocode.clades.confident_clades — без адаптации, как договорились.
    """
    newick = treefile.read_text(encoding="utf-8").strip()
    if not newick:
        return []
    tree = TreeResult(method="iqtree", newick=newick)
    return clades_mod.confident_clades(tree, muts=[], ufboot_min=ufboot_min, alrt_min=alrt_min)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--mrbayes-dir", default="mrbayes",
                    help="папка с <группа>.nex.con.tre (+ .names.tsv); '' — отключить источник")
    ap.add_argument("--iqtree-dir", default="",
                    help="папка вида trees/<группа>/<группа>.treefile; по умолчанию отключено")
    ap.add_argument("--posterior-min", type=float, default=0.95)
    ap.add_argument("--ufboot-min", type=float, default=95.0)
    ap.add_argument("--alrt-min", type=float, default=80.0)
    ap.add_argument("--out", default="groups/report.json")
    args = ap.parse_args(argv)

    report: dict[str, dict] = {}

    if args.mrbayes_dir:
        mb_dir = Path(args.mrbayes_dir)
        if mb_dir.is_dir():
            for con_tre in sorted(mb_dir.glob(f"*{CON_TRE_SUFFIX}")):
                key = con_tre.name[: -len(CON_TRE_SUFFIX)]
                names_tsv = mb_dir / f"{key}.names.tsv"
                cl = clades_from_mrbayes(con_tre, names_tsv, args.posterior_min)
                report.setdefault(key, {})["mrbayes"] = {
                    "threshold": {"posterior_min": args.posterior_min}, "clades": cl}
                print(f"[{key}] mrbayes: {len(cl)} уверенных клад (posterior≥{args.posterior_min})")
        else:
            print(f"--mrbayes-dir не найдена: {mb_dir}", file=sys.stderr)

    if args.iqtree_dir:
        iq_dir = Path(args.iqtree_dir)
        if iq_dir.is_dir():
            for sub in sorted(p for p in iq_dir.iterdir() if p.is_dir()):
                treefile = sub / f"{sub.name}.treefile"
                if not treefile.is_file():
                    continue
                cl = clades_from_iqtree(treefile, args.ufboot_min, args.alrt_min)
                report.setdefault(sub.name, {})["iqtree"] = {
                    "threshold": {"ufboot_min": args.ufboot_min, "alrt_min": args.alrt_min},
                    "clades": cl}
                print(f"[{sub.name}] iqtree: {len(cl)} уверенных клад "
                      f"(UFBoot≥{args.ufboot_min}, aLRT≥{args.alrt_min})")
        else:
            print(f"--iqtree-dir не найдена: {iq_dir}", file=sys.stderr)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")

    n_clades = sum(len(v.get("mrbayes", {}).get("clades", [])) + len(v.get("iqtree", {}).get("clades", []))
                   for v in report.values())
    print(f"\n{out_path}: {len(report)} групп, {n_clades} уверенных клад суммарно")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
