#!/usr/bin/env python3
"""HTML/SVG-визуализация дерева — ОТДЕЛЬНО от тестового рантайма.

Намеренно вынесено из test_pipeline.py / test_fixtures.py: те тестируют
pipeline через консоль (subprocess + assert на report.json/.nex.con.tre) и не
импортируют этот файл вообще — так что консольные тесты работают, даже если
здесь что-то сломано или Bio.Phylo/HTML-генерация недоступны по любой причине.

Использование как CLI (после того как build_trees_mrbayes.py и
clade_search.py уже отработали):

    python visualize_tree.py <group_key> <mrbayes_dir> [--report report.json] [--out out.html]

Пример:
    python visualize_tree.py IGHV3-23_01_IGHJ3_01 ../mrbayes \\
        --report ../clades/report.json --out ../mrbayes/IGHV3-23_01_IGHJ3_01.tree.html

Или программно:
    from visualize_tree import parse_nexus_tree, generate_tree_html
"""
from __future__ import annotations

import argparse
import json
import sys
from io import StringIO
from pathlib import Path


def parse_nexus_tree(nex_path: Path, names_map_path: Path | None = None) -> tuple[str | None, dict]:
    """Парсить .nex.con.tre (NEXUS консенсус) → (newick, {clade_sig: {posterior}})."""
    from Bio import Phylo

    if not nex_path.is_file():
        return None, {}

    rev = {}
    if names_map_path and names_map_path.is_file():
        for line in names_map_path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            safe, orig = line.split("\t", 1)
            rev[safe] = orig

    try:
        tree = next(Phylo.parse(str(nex_path), "nexus"), None)
        if tree is None:
            return None, {}

        if rev:
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
    except Exception as e:
        print(f"  ⚠ Ошибка парсинга {nex_path}: {e}", file=sys.stderr)
        return None, {}


def generate_tree_html(group_key: str, newick: str, clades: list[dict], out_path: Path) -> bool:
    """Генерировать HTML с филограммой (rectangular phylogram) на основе newick и клад.

    Возвращает True при успехе.
    """
    from Bio import Phylo

    try:
        tree = Phylo.read(StringIO(newick), "newick")
    except Exception as e:
        print(f"  ⚠ Ошибка парсинга newick {group_key}: {e}", file=sys.stderr)
        return False

    tips = tree.get_terminals()
    n_tips = len(tips)
    if n_tips == 0:
        return False

    leaf_order = []

    def dfs(clade):
        if clade.is_terminal():
            leaf_order.append(clade.name)
        else:
            for child in clade.clades:
                dfs(child)

    dfs(tree.root)

    def max_depth_of(clade, depth=0):
        if clade.is_terminal():
            return depth
        return max(max_depth_of(c, depth + (c.branch_length or 0)) for c in clade.clades) \
               if clade.clades else depth

    max_d = max_depth_of(tree.root) or 1

    left_margin = 40
    plot_w = 520
    scale = plot_w / max_d
    row_h = 56
    top = 60

    y_coords = {t: top + i * row_h for i, t in enumerate(leaf_order)}

    def x_coord(depth: float) -> float:
        return left_margin + depth * scale

    def build_coords(clade, depth=0):
        x = x_coord(depth)
        if clade.is_terminal():
            return {"name": clade.name, "x": x, "y": y_coords.get(clade.name, 0),
                    "type": "leaf", "confidence": clade.confidence}
        children = [build_coords(c, depth + (c.branch_length or 0)) for c in clade.clades]
        avg_y = sum(c["y"] for c in children) / len(children) if children else 0
        return {"x": x, "y": avg_y, "type": "internal", "children": children,
                "confidence": clade.confidence}

    tree_coords = build_coords(tree.root)

    confident_leaf_pairs = {tuple(sorted(c["leaves"]))
                           for c in clades if c.get("posterior", 0) >= 0.95}

    svg_lines = [
        f'<svg viewBox="0 0 760 {top + len(leaf_order) * row_h + 60}" '
        f'role="img" aria-label="Дерево {group_key}">',
    ]

    def get_leaves(node):
        if node["type"] == "leaf":
            return [node["name"]]
        return sum((get_leaves(c) for c in node.get("children", [])), [])

    def render_branches(node, parent_x=None, parent_y=None, parent_confident=False):
        if node["type"] == "leaf":
            x, y = node["x"], node["y"]
            if parent_x is not None:
                is_confident = any(node["name"] in pair for pair in confident_leaf_pairs)
                css_class = "branch confident" if is_confident else "branch"
                svg_lines.append(f'<line class="{css_class}" x1="{parent_x}" y1="{parent_y}" '
                               f'x2="{x}" y2="{y}"></line>')
            css_class = "leaf-dot confident" if parent_confident else "leaf-dot"
            svg_lines.append(f'<circle class="{css_class}" cx="{x}" cy="{y}" r="3.5"></circle>')
            svg_lines.append(f'<circle class="hit" cx="{x}" cy="{y}" r="10">'
                           f'<title>{node["name"]}</title></circle>')
            label_short = node["name"][:14] + ("…" if len(node["name"]) > 14 else "")
            svg_lines.append(f'<text class="leaf-label" x="{x+8}" y="{y-3}">{label_short}</text>')
        else:
            x, y = node["x"], node["y"]
            if parent_x is not None:
                svg_lines.append(f'<line class="branch" x1="{parent_x}" y1="{parent_y}" '
                               f'x2="{x}" y2="{parent_y}"></line>')
                svg_lines.append(f'<line class="branch" x1="{x}" y1="{parent_y}" '
                               f'x2="{x}" y2="{y}"></line>')
            if node["children"]:
                y_min = min(c["y"] for c in node["children"])
                y_max = max(c["y"] for c in node["children"])
                svg_lines.append(f'<line class="branch" x1="{x}" y1="{y_min}" '
                               f'x2="{x}" y2="{y_max}"></line>')
            if node["confidence"] is not None:
                pp = float(node["confidence"])
                pp = pp / 100 if pp > 1 else pp
                svg_lines.append(f'<text class="node-label" x="{x-5}" y="{y-5}" text-anchor="end">'
                               f'PP {round(pp, 2)}</text>')
            for child in node.get("children", []):
                child_is_conf = tuple(sorted(get_leaves(child))) in confident_leaf_pairs
                render_branches(child, x, y, child_is_conf)

    render_branches(tree_coords)
    svg_lines.append('</svg>')

    n_confident = len([c for c in clades if c.get("posterior", 0) >= 0.95])
    html = f"""<title>Дерево {group_key}</title>
<style>
  .viz-root {{
    color-scheme: light;
    --surface-1:      #fcfcfb;
    --page:           #f9f9f7;
    --text-primary:   #0b0b0b;
    --text-secondary: #52514e;
    --text-muted:     #898781;
    --gridline:       #e1e0d9;
    --branch:         #52514e;
    --good:           #0ca30c;
    --border:         rgba(11,11,11,0.10);
  }}
  @media (prefers-color-scheme: dark) {{
    :root:where(:not([data-theme="light"])) .viz-root {{
      color-scheme: dark;
      --surface-1:      #1a1a19;
      --page:           #0d0d0d;
      --text-primary:   #ffffff;
      --text-secondary: #c3c2b7;
      --text-muted:     #898781;
      --gridline:       #2c2c2a;
      --branch:         #c3c2b7;
      --good:           #0ca30c;
      --border:         rgba(255,255,255,0.10);
    }}
  }}
  * {{ box-sizing: border-box; }}
  body {{ margin: 0; }}
  .viz-root {{
    background: var(--page); color: var(--text-primary);
    font-family: system-ui, -apple-system, "Segoe UI", sans-serif;
    padding: 28px 20px;
  }}
  .wrap {{ max-width: 860px; margin: 0 auto; }}
  .card {{
    background: var(--surface-1); border: 1px solid var(--border);
    border-radius: 14px; padding: 22px 24px; margin-bottom: 16px;
  }}
  h1 {{ font-size: 18px; margin: 0 0 4px; }}
  .meta {{
    display: flex; flex-wrap: wrap; gap: 8px 18px;
    font-size: 12.5px; color: var(--text-muted);
    margin-top: 10px; padding-top: 10px; border-top: 1px solid var(--gridline);
  }}
  .meta b {{ color: var(--text-secondary); font-weight: 600; }}
  svg {{ display: block; width: 100%; height: auto; }}
  .branch {{ stroke: var(--branch); stroke-width: 2; fill: none; stroke-linecap: round; }}
  .branch.confident {{ stroke: var(--good); stroke-width: 2.5; }}
  .leaf-dot {{ fill: var(--text-primary); }}
  .leaf-dot.confident {{ fill: var(--good); }}
  .hit {{ fill: transparent; cursor: default; }}
  .leaf-label {{ font-size: 12px; fill: var(--text-primary); dominant-baseline: middle; }}
  .node-label {{ font-size: 11px; fill: var(--good); font-weight: 600; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 12.5px; margin-top: 6px; }}
  th, td {{ text-align: left; padding: 6px 10px; border-bottom: 1px solid var(--gridline); }}
  th {{ color: var(--text-muted); font-weight: 600; font-size: 11.5px; }}
  td.n {{ text-align: right; font-variant-numeric: tabular-nums; }}
  .pill {{
    display: inline-block; padding: 1px 8px; border-radius: 8px;
    background: color-mix(in srgb, var(--good) 16%, transparent);
    color: var(--good); font-weight: 600; font-size: 11.5px;
  }}
</style>
<div class="viz-root">
  <div class="wrap">
    <div class="card">
      <h1>Дерево {group_key}</h1>
      <div class="meta">
        <span><b>таксонов</b> {n_tips}</span>
        <span><b>уверенных клад</b> {n_confident}</span>
      </div>
    </div>
    <div class="card">
      {"".join(svg_lines)}
    </div>
"""
    if n_confident:
        html += """    <div class="card">
      <h1 style="font-size:15px;">Уверенные клады (posterior ≥ 0.95)</h1>
      <table>
        <tr><th>листья</th><th class="n">posterior</th></tr>
"""
        for c in clades:
            if c.get("posterior", 0) >= 0.95:
                leaves_str = ", ".join(c["leaves"][:2]) + ("..." if len(c["leaves"]) > 2 else "")
                html += (f'        <tr><td><span class="pill">✓</span> {leaves_str}</td>'
                       f'<td class="n">{c.get("posterior", "—")}</td></tr>\n')
        html += "      </table>\n    </div>\n"

    html += "  </div>\n</div>"
    out_path.write_text(html, encoding="utf-8")
    return True


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Визуализация дерева группы как HTML.")
    ap.add_argument("group_key", help="ключ группы (напр. IGHV3-23_01_IGHJ3_01)")
    ap.add_argument("mrbayes_dir", type=Path, help="папка с <группа>.nex.con.tre (+.names.tsv)")
    ap.add_argument("--report", type=Path, default=None,
                    help="clades/report.json (для подсветки уверенных клад)")
    ap.add_argument("--out", type=Path, default=None, help="путь выходного .html")
    args = ap.parse_args(argv)

    con_tre = args.mrbayes_dir / f"{args.group_key}.nex.con.tre"
    names_tsv = args.mrbayes_dir / f"{args.group_key}.names.tsv"
    newick, _ = parse_nexus_tree(con_tre, names_tsv)
    if not newick:
        print(f"Не удалось получить дерево из {con_tre}", file=sys.stderr)
        return 1

    clades = []
    if args.report and args.report.is_file():
        report = json.loads(args.report.read_text(encoding="utf-8"))
        clades = report.get(args.group_key, {}).get("mrbayes", {}).get("clades", [])

    out_path = args.out or (args.mrbayes_dir / f"{args.group_key}.tree.html")
    ok = generate_tree_html(args.group_key, newick, clades, out_path)
    if ok:
        print(f"✓ {out_path}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
