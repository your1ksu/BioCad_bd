#!/usr/bin/env python3
"""ШАГ «Никиты» (уверенные клады) — nexus → report.json.

Обёртка над продакшн-движком ``biocode`` (см. ../biocode/, вендорено без
изменений из BIOCAD.bigchallenges@main), адаптированная под файловый контракт
остальных участников конвейера (см. ../BioCad_repo.md, строка «уверенные
клады»: критерий UFBoot≥95 и aLRT≥80, у Байеса — posterior≥0.95).

Два независимых источника NEXUS-деревьев поддержаны (можно использовать один
или оба сразу — тогда клады группы объединяются в один report.json):

  --mrbayes-dir DIR (по умолчанию 'mrbayes')
      <группа>.nex.con.tre из шага ../04b_build_trees_mrbayes/build_trees_mrbayes.py — консенсус
      MrBayes с апостериорной поддержкой (posterior) на внутренних узлах.
      Критерий: posterior ≥ --posterior-min (по умолчанию 0.95).
      Требуется <группа>.names.tsv рядом (тот же шаг его пишет), чтобы вернуть
      исходные id таксонов вместо MrBayes-совместимых T0001...

      Два исправления против наивного обхода дерева (найдены и проверены
      тестами на tests/fixtures/, см. tests/test_fixtures.py):
      1) MrBayes пишет неукоренённое дерево как rooted-newick с базовой
         политомией — часть клады может остаться «голыми» листьями прямо на
         корне вместо отдельного узла (см. _root_complement_supports).
      2) Обе стороны одного и того же ребра дерева иногда попадают в supports
         как два разных «клады» — оставляем только содержательную меньшую
         сторону (см. _drop_complement_duplicates).

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
from io import StringIO
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from biocode import clades as clades_mod
from biocode.model import TreeResult
from biocode.trees import mrbayes as mrbayes_mod

CON_TRE_SUFFIX = ".nex.con.tre"


def _root_complement_supports(newick: str, supports: dict[str, dict]) -> dict[str, dict]:
    """Достроить supports той стороной корневого разбиения, которую MrBayes не
    обернул в отдельный узел.

    MrBayes пишет консенсус НЕУКОРЕНЁННОГО дерева как rooted-newick с базовой
    политомией (напр. ``(seq1,seq2,(seq3,seq4)1.00)``): 2 таксона висят прямо
    на корне по отдельности, а их сестринская пара обёрнута в узел с posterior.
    Для 4 таксонов это ЕДИНСТВЕННОЕ возможное внутреннее ребро дерева — то есть
    {seq1,seq2} обладает ТОЙ ЖЕ posterior-поддержкой, что и {seq3,seq4}, просто
    не материализована как узел, и наивный обход get_nonterminals() её теряет.

    Правило (безопасно и однозначно только в этом частном случае): если среди
    прямых потомков КОРНЯ ровно один — свёрнутая клада, а остальные ≥2 —
    голые листья, то эти листья вместе — комплементарная клада с той же
    posterior (это то же самое ребро дерева, две стороны одного разбиения).
    Если на корне ≥2 свёрнутых клад — топология между ними не разрешена,
    объединять голые листья не тождественно правильно, поэтому не трогаем.
    """
    if not newick:
        return supports
    try:
        from Bio import Phylo
        tree = Phylo.read(StringIO(newick), "newick")
    except Exception:
        return supports

    root_children = list(tree.root.clades)
    internal_children = [c for c in root_children if not c.is_terminal()]
    leaf_children = [c for c in root_children if c.is_terminal()]
    if len(internal_children) != 1 or len(leaf_children) < 2:
        return supports

    wrapped = internal_children[0]
    wrapped_leaves = frozenset(t.name for t in wrapped.get_terminals())
    wrapped_sig = "|".join(sorted(wrapped_leaves))
    if wrapped_sig not in supports:
        return supports

    complement_leaves = sorted(t.name for t in leaf_children)
    complement_sig = "|".join(complement_leaves)
    if complement_sig in supports:
        return supports

    out = dict(supports)
    out[complement_sig] = dict(supports[wrapped_sig])
    return out


def _drop_complement_duplicates(supports: dict[str, dict], all_leaves: frozenset) -> dict[str, dict]:
    """Если клада C и её дополнение (all_leaves \\ C) обе есть в supports — это
    ДВЕ СТОРОНЫ ОДНОГО РЕБРА дерева. Когда стороны разного размера, большая —
    это просто «всё, кроме меньшей клады» и не несёт отдельного биологического
    смысла (пример: клада из 2 листьев vs комплемент из 4 разнородных) — такую
    большую сторону убираем, оставляя только содержательную меньшую.

    Когда стороны РАВНОГО размера (напр. 2+2 при 4 таксонах) — это два
    ОДИНАКОВО специфичных, независимо содержательных разбиения (напр. {1,2} и
    {3,4} — обе реальные пары), обе сохраняются.
    """
    sig_to_leaves = {sig: frozenset(sig.split("|")) for sig in supports}
    drop: set[str] = set()
    sigs = list(supports)
    for i, sig_a in enumerate(sigs):
        if sig_a in drop:
            continue
        leaves_a = sig_to_leaves[sig_a]
        complement = all_leaves - leaves_a
        if len(complement) < 2 or len(leaves_a) == len(complement):
            continue
        for sig_b in sigs[i + 1:]:
            if sig_b in drop or sig_to_leaves[sig_b] != complement:
                continue
            drop.add(sig_a if len(leaves_a) > len(complement) else sig_b)
            break
    return {sig: sup for sig, sup in supports.items() if sig not in drop}


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
    newick, supports = mrbayes_mod.parse_con_tre(con_tre, rev)
    supports = _root_complement_supports(newick, supports)

    if newick:
        from Bio import Phylo
        all_leaves = frozenset(t.name for t in Phylo.read(StringIO(newick), "newick").get_terminals())
        supports = _drop_complement_duplicates(supports, all_leaves)

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


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    key = None
    seq: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if key is not None:
                records[key] = "".join(seq)
            key = line[1:].split()[0]
            seq = []
        else:
            seq.append(line)
    if key is not None:
        records[key] = "".join(seq)
    return records


def _find_aligned_fasta(group_key: str, aligned_dir: Path) -> Path | None:
    candidates = [
        aligned_dir / f"{group_key}.fasta",
        aligned_dir / f"{group_key}.fa",
        aligned_dir / f"{group_key}_aligned.fasta",
        aligned_dir / f"{group_key}_aligned.fa",
    ]
    for p in candidates:
        if p.is_file():
            return p
    return None


def write_clade_fastas(report: dict[str, dict],
                       aligned_dir: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for group_key, sources in report.items():
        aligned = _find_aligned_fasta(group_key, aligned_dir)
        if aligned is None:
            print(f"  {group_key}: aligned FASTA не найден")
            continue
        all_seqs = read_fasta(aligned)
        for src_key in ("mrbayes", "iqtree"):
            for ci, clade in enumerate(sources.get(src_key, {}).get("clades", [])):
                leaves = clade.get("leaves", [])
                leaf_seqs = {k: v for k, v in all_seqs.items() if k in leaves}
                if not leaf_seqs:
                    continue
                fasta_path = out_dir / f"{group_key}__{src_key}_c{ci:03d}.fa"
                lines = [f">{k}\n{v}" for k, v in leaf_seqs.items()]
                fasta_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        n_clades = (len(sources.get("mrbayes", {}).get("clades", []))
                    + len(sources.get("iqtree", {}).get("clades", [])))
        print(f"  {group_key}: {n_clades} clade FASTA files")


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
    ap.add_argument("--aligned-dir", default="",
                    help="папка с выравненными FASTA (группа.fasta); если указан, извлекаются FASTA клад")
    ap.add_argument("--clades-fasta-dir", default="",
                    help="куда сохранять FASTA уверенных клад (требует --aligned-dir)")
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

    if args.aligned_dir and args.clades_fasta_dir:
        write_clade_fastas(report, Path(args.aligned_dir), Path(args.clades_fasta_dir))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
