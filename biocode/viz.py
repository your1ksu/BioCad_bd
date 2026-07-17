"""Визуализация деревьев и FR/CDR (блок B12).

Рендер филограммы группы (matplotlib, headless):
  • листья раскрашены по изотипу (IGHM/IGHD/IGHG/IGHA…), germline-outgroup помечен;
  • ветви окрашены по числу мутаций (тепловая интенсивность);
  • уверенные клады подсвечены и подписаны (UFBoot/aLRT, при Bayes — posterior);
плюс цветная полоса FR/CDR по колонкам выравнивания.
"""
from __future__ import annotations

from collections import Counter
from io import StringIO
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from .model import Mutation, SequenceRecord, TreeResult

ISOTYPE_COLOR = {"IGHM": "#7c93b3", "IGHD": "#a9b8cf", "IGHG": "#d1553e",
                 "IGHA": "#e8a23c", "IGHE": "#8b56ce", None: "#c9c9c9"}
REGION_COLOR = {"FR1": "#e8edf3", "FR2": "#dbe3ec", "FR3": "#cfd9e6", "FR4": "#c3cfe0",
                "CDR1": "#f4c542", "CDR2": "#f08a3c", "CDR3": "#d94f3d", "other": "#ffffff"}
REGION_ORDER = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]


def _node_name(clade) -> str | None:
    if clade.is_terminal():
        return clade.name
    return clade.name.split("/")[0] if clade.name else None


def _support(clade) -> dict:
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
    return {"ufboot": nums[0]} if nums else {}


def plot_tree(tree: TreeResult, records: dict[str, SequenceRecord],
              muts: list[Mutation], clades: list[dict], out_path: str | Path,
              *, ufboot_min: float = 95.0, alrt_min: float = 80.0) -> Path:
    """Филограмма группы с раскраской по изотипу, мутациям и уверенным кладам."""
    from Bio import Phylo

    tree_obj = Phylo.read(StringIO(tree.newick), "newick")
    tree_obj.ladderize()
    terminals = tree_obj.get_terminals()
    y_of = {id(t): i for i, t in enumerate(terminals)}

    # x (глубина по длинам ветвей) и y (внутренние = среднее детей)
    x_of: dict[int, float] = {}

    def assign_x(clade, x0):
        x = x0 + (clade.branch_length or 0.0)
        x_of[id(clade)] = x
        for c in clade.clades:
            assign_x(c, x)
    assign_x(tree_obj.root, 0.0)

    def assign_y(clade):
        if clade.is_terminal():
            return y_of[id(clade)]
        ys = [assign_y(c) for c in clade.clades]
        y = sum(ys) / len(ys)
        y_of[id(clade)] = y
        return y
    assign_y(tree_obj.root)

    parent_of: dict[int, object] = {}
    for clade in tree_obj.find_clades():
        for c in clade.clades:
            parent_of[id(c)] = clade

    # мутаций на входящую ветвь узла (для окраски ветвей)
    mut_by_child: Counter = Counter(m.branch.split("→")[-1] for m in muts)
    max_mut = max(mut_by_child.values(), default=1)

    fig, ax = plt.subplots(figsize=(11, max(4, 0.32 * len(terminals) + 2)),
                           facecolor="white")

    # ветви
    for clade in tree_obj.find_clades():
        cx, cy = x_of[id(clade)], y_of[id(clade)]
        if clade is not tree_obj.root:
            parent = parent_of.get(id(clade))
            if parent is not None:
                px = x_of[id(parent)]
                nk = _node_name(clade)
                load = mut_by_child.get(nk or clade.name, 0)
                shade = 0.2 + 0.6 * (load / max_mut)
                col = (0.85 - shade * 0.6, 0.2, 0.2) if load else (0.6, 0.6, 0.6)
                lw = 0.8 + 2.2 * (load / max_mut)
                ax.plot([px, cx], [cy, cy], color=col, lw=lw, solid_capstyle="round", zorder=2)
        if not clade.is_terminal():
            ys = [y_of[id(c)] for c in clade.clades]
            ax.plot([cx, cx], [min(ys), max(ys)], color="#888", lw=1, zorder=1)

    # уверенные клады — подсветка
    conf_nodes = {c["clade"] for c in clades}
    for clade in tree_obj.get_nonterminals():
        nk = _node_name(clade)
        if nk in conf_nodes:
            leaves = clade.get_terminals()
            ys = [y_of[id(t)] for t in leaves]
            ax.axhspan(min(ys) - 0.4, max(ys) + 0.4, xmin=0, xmax=1,
                       color="#42a463", alpha=0.06, zorder=0)
            sup = _support(clade)
            ci = next((c for c in clades if c["clade"] == nk), {})
            lbl = f"{nk}: UF{sup.get('ufboot',0):.0f}"
            if ci.get("posterior") is not None:
                lbl += f"/pp{ci['posterior']:.2f}"
            ax.text(x_of[id(clade)], max(ys) + 0.5, lbl, fontsize=6.5,
                    color="#2f7d4e", va="bottom")

    # листья
    xmax = max(x_of.values()) if x_of else 1
    for t in terminals:
        rec = records.get(t.name)
        iso = rec.isotype if rec else None
        is_germ = t.name.startswith("GERMLINE_")
        color = "#111111" if is_germ else ISOTYPE_COLOR.get(iso, ISOTYPE_COLOR[None])
        marker = "s" if is_germ else "o"
        ax.scatter([x_of[id(t)]], [y_of[id(t)]], s=26, color=color,
                   edgecolors="#333", linewidth=0.4, marker=marker, zorder=3)
        name = "germline" if is_germ else t.name[:18]
        ax.text(x_of[id(t)] + xmax * 0.01, y_of[id(t)], name, fontsize=6,
                va="center", color="#222")

    isos = [i for i in ISOTYPE_COLOR if i and any(
        (records.get(t.name).isotype if records.get(t.name) else None) == i for t in terminals)]
    handles = [Patch(color=ISOTYPE_COLOR[i], label=i) for i in isos]
    handles += [Patch(color="#111111", label="germline"),
                Patch(color="#d1553e", label="ветвь: мутации (интенсивность)"),
                Patch(color="#42a463", alpha=0.2, label="уверенная клада")]
    ax.legend(handles=handles, loc="lower right", fontsize=7, framealpha=0.9)

    ax.set_title(f"{tree.method.upper()} дерево · model={tree.model} · "
                 f"{len(terminals)} листьев · уверенных клад: {len(clades)}",
                 fontsize=11, fontweight="bold")
    ax.set_xlabel("эволюционная дистанция (замен/сайт)", fontsize=9)
    ax.set_yticks([])
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.set_xlim(-xmax * 0.02, xmax * 1.25)
    ax.set_ylim(-1, len(terminals))
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


def plot_region_track(track: list[dict], out_path: str | Path) -> Path:
    """Цветная полоса FR/CDR по колонкам выравнивания."""
    fig, ax = plt.subplots(figsize=(11, 1.1), facecolor="white")
    for t in track:
        ax.axvspan(t["column"], t["column"] + 1,
                   color=REGION_COLOR.get(t["region"], "#fff"))
    ax.set_xlim(0, len(track))
    ax.set_yticks([])
    ax.set_xlabel("колонка выравнивания", fontsize=9)
    ax.set_title("FR/CDR-разметка", fontsize=10, fontweight="bold")
    handles = [Patch(color=REGION_COLOR[r], label=r) for r in REGION_ORDER]
    ax.legend(handles=handles, loc="upper center", ncol=7, fontsize=7,
              bbox_to_anchor=(0.5, -0.4), framealpha=0)
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path
