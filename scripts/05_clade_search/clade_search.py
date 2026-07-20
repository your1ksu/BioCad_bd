#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from io import StringIO
from pathlib import Path

CON_TRE_SUFFIX = ".nex.con.tre"


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


def _clade_tree_metrics(tree, leaves: frozenset):
    """Метрики клады по дереву: (depth, ancestor_to_leaves).

    depth — сколько ВНУТРЕННИХ узлов внутри клады, не считая сам узел-предок
    клады. Пара листьев от общего предка → depth=0; каждый доп. уровень
    ветвления внутри клады → +1. Считается по РЕАЛЬНОМУ дереву, а не как
    (size-1): majority-rule консенсус MrBayes может давать политомии.

    ancestor_to_leaves — насколько предок клады далеко от своих листьев, в
    длинах ветвей (замен на сайт): {"max","mean","min"} по всем листьям клады.

    Клада может быть достроена из комплемента корня (см.
    _root_complement_supports): тогда отдельного узла в дереве нет, а общим
    предком листьев выступает корень — обрабатывается отдельной веткой.
    """
    if tree is None:
        return None, None
    target = None
    for clade in tree.get_nonterminals():
        if frozenset(t.name for t in clade.get_terminals()) == leaves:
            target = clade
            break

    if target is not None:
        depth = len(target.get_nonterminals()) - 1        # без самого узла-предка
        # Bio.Phylo depths() у под-клады стартует с её собственной ветви (ведущей
        # К предку) — вычитаем её, чтобы получить расстояние ОТ предка вниз.
        stem = target.branch_length or 0.0
        dists = [d - stem for cl, d in target.depths().items() if cl.is_terminal()]
    else:
        # комплемент корня: листья висят прямо на корне, их предок — корень
        root_leaf_children = [c for c in tree.root.clades
                              if c.is_terminal() and c.name in leaves]
        if len(root_leaf_children) != len(leaves):
            return None, None
        depth = 0
        dists = [c.branch_length or 0.0 for c in root_leaf_children]

    if not dists:
        return depth, None
    return depth, {
        "max": round(max(dists), 6),
        "mean": round(sum(dists) / len(dists), 6),
        "min": round(min(dists), 6),
    }


def confident_clades(newick: str, *, ufboot_min: float = 95.0,
                     alrt_min: float = 80.0) -> list[dict]:
    from Bio import Phylo

    if not newick:
        return []
    phylo = Phylo.read(StringIO(newick), "newick")

    clades: list[dict] = []
    for clade in phylo.get_nonterminals():
        sup = _support(clade)
        uf, al = sup.get("ufboot", 0.0), sup.get("alrt", 0.0)
        if not (uf >= ufboot_min and al >= alrt_min):
            continue
        leaves = [t.name for t in clade.get_terminals()]
        if len(leaves) < 2:
            continue
        depth, anc = _clade_tree_metrics(phylo, frozenset(leaves))
        clades.append({
            "clade": _node_name(clade),
            "size": len(leaves),
            "leaves": leaves,
            "depth": depth,
            "ancestor_to_leaves": anc,
            "ufboot": uf,
            "alrt": al,
            "defining_mutations": 0,
            "defining_cdr": 0,
            "isotypes": {},
            "days": [],
            "confident_both_models": False,
        })
    clades.sort(key=lambda c: (c["ufboot"], c["size"]), reverse=True)
    return clades


def _root_complement_supports(newick: str, supports: dict[str, dict]) -> dict[str, dict]:
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


def parse_con_tre(path: str | Path, rev: dict[str, str]) -> tuple[str, dict[str, dict]]:
    from Bio import Phylo

    path = Path(path)
    if not path.is_file():
        return "", {}
    tree = next(Phylo.parse(str(path), "nexus"), None)
    if tree is None:
        return "", {}
    for tip in tree.get_terminals():
        tip.name = rev.get(tip.name, tip.name)

    all_leaves = {t.name for t in tree.get_terminals()}
    supports: dict[str, dict] = {}
    for clade in tree.get_nonterminals():
        leaves = frozenset(t.name for t in clade.get_terminals())
        if len(leaves) < 2 or len(leaves) >= len(all_leaves):
            continue
        if clade.confidence is not None:
            pp = float(clade.confidence)
            pp = pp / 100.0 if pp > 1.0 else pp
            supports["|".join(sorted(leaves))] = {"posterior": round(pp, 4)}
    out = StringIO()
    Phylo.write(tree, out, "newick")
    return out.getvalue().strip(), supports


def _load_names(names_tsv: Path) -> dict[str, str]:
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
    rev = _load_names(names_tsv)
    newick, supports = parse_con_tre(con_tre, rev)
    supports = _root_complement_supports(newick, supports)

    tree = None
    if newick:
        from Bio import Phylo
        tree = Phylo.read(StringIO(newick), "newick")
        all_leaves = frozenset(t.name for t in tree.get_terminals())
        supports = _drop_complement_duplicates(supports, all_leaves)

    entries = []
    for sig, sup in supports.items():
        pp = sup.get("posterior")
        if pp is None or pp < posterior_min:
            continue
        leaves = sig.split("|")
        depth, anc = _clade_tree_metrics(tree, frozenset(leaves))
        entries.append({
            "clade": sig, "size": len(leaves), "leaves": leaves,
            "depth": depth, "ancestor_to_leaves": anc,
            "ufboot": None, "alrt": None, "posterior": pp,
            "defining_mutations": None, "defining_cdr": None,
            "isotypes": {}, "days": [], "confident_both_models": False,
        })
    entries.sort(key=lambda c: (c["posterior"], c["size"]), reverse=True)
    return entries


def clades_from_iqtree(treefile: Path, ufboot_min: float, alrt_min: float) -> list[dict]:
    newick = treefile.read_text(encoding="utf-8").strip()
    if not newick:
        return []
    return confident_clades(newick, ufboot_min=ufboot_min, alrt_min=alrt_min)


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
                clean = {k: v.replace("-", "").replace(".", "") for k, v in leaf_seqs.items()}
                lines = [f">{k}\n{v}" for k, v in clean.items()]
                fasta_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        n_clades = (len(sources.get("mrbayes", {}).get("clades", []))
                    + len(sources.get("iqtree", {}).get("clades", [])))
        print(f"  {group_key}: {n_clades} clade FASTA files")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Уверенные клады — ML + Bayes")
    ap.add_argument("--mrbayes-dir", default="mrbayes",
                    help="папка с <группа>.nex.con.tre (+ .names.tsv)")
    ap.add_argument("--iqtree-dir", default="",
                    help="папка вида trees/<группа>/<группа>.treefile")
    ap.add_argument("--posterior-min", type=float, default=0.95)
    ap.add_argument("--ufboot-min", type=float, default=95.0)
    ap.add_argument("--alrt-min", type=float, default=80.0)
    ap.add_argument("--out", default="groups/report.json")
    ap.add_argument("--aligned-dir", default="",
                    help="папка с выравненными FASTA")
    ap.add_argument("--clades-fasta-dir", default="",
                    help="куда сохранять FASTA уверенных клад")
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

    # кросс-модельное подтверждение: клада надёжна в обеих моделях, если её набор
    # листьев присутствует и в mrbayes, и в iqtree одной группы
    for sources in report.values():
        mb_sets = {frozenset(c["leaves"]) for c in sources.get("mrbayes", {}).get("clades", [])}
        iq_sets = {frozenset(c["leaves"]) for c in sources.get("iqtree", {}).get("clades", [])}
        for c in sources.get("mrbayes", {}).get("clades", []):
            c["confident_both_models"] = frozenset(c["leaves"]) in iq_sets
        for c in sources.get("iqtree", {}).get("clades", []):
            c["confident_both_models"] = frozenset(c["leaves"]) in mb_sets

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
