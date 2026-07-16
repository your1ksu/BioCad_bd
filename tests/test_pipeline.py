#!/usr/bin/env python3
"""Полный тестовый цикл: запуск скриптов + генерация визуализаций.

Этот скрипт:
1. Загружает несколько тестовых групп из BIOCAD.bigchallenges (реальные данные)
2. Запускает mrbayes/run_mrbayes.py
3. Запускает groups/confident_clades_report.py
4. Парсит результаты (NEXUS, JSON)
5. Генерирует HTML-визуализации филограмм для каждой группы
6. Выводит статистику и пути к результатам
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from io import StringIO

# Требует: biopython (для парсинга NEXUS)
try:
    from Bio import Phylo
except ImportError:
    print("Ошибка: требуется pybiopython. Установка: pip install biopython", file=sys.stderr)
    sys.exit(1)


def download_test_groups(n_groups: int = 3) -> dict[str, Path]:
    """Загрузить несколько тестовых групп из BIOCAD.bigchallenges.

    Возвращает {group_key: fasta_path} для групп разного размера.
    """
    repo_path = Path(__file__).parent / "BIOCAD.bigchallenges"
    if not repo_path.is_dir():
        print(f"Ошибка: репозиторий не найден: {repo_path}", file=sys.stderr)
        return {}

    # Найти доступные части (part_1, part_2, part_3, part_4)
    anotherp = repo_path / "anotherpipeline"
    available_files = []
    for part in ["part_1", "part_2", "part_3", "part_4"]:
        part_dir = anotherp / part
        if part_dir.is_dir():
            for fasta in sorted(part_dir.glob("*_aligned.fasta"))[:20]:  # берём первые 20 файлов
                group_key = fasta.stem.replace("_aligned", "")
                available_files.append((group_key, fasta))

    in_dir = Path(__file__).parent / "aligned_sequences"
    in_dir.mkdir(exist_ok=True)
    result = {}

    for group_key, src in available_files[:n_groups]:
        dst = in_dir / f"{group_key}_aligned.fasta"
        dst.write_bytes(src.read_bytes())
        result[group_key] = dst
        n_seqs = src.read_text().count(">")
        print(f"  ✓ {group_key}: {n_seqs} таксонов")

    return result


def parse_nexus_tree(nex_path: Path, names_map_path: Path | None = None) -> tuple[str | None, dict]:
    """Парсить .nex.con.tre (NEXUS консенсус) → (newick, {clade_sig: {posterior}})."""
    if not nex_path.is_file():
        return None, {}

    # Если есть .names.tsv, загрузить обратный маппинг (safe_id → original_id)
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

        # Восстановить исходные id таксонов
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


def generate_tree_html(group_key: str, newick: str, clades: list[dict],
                       mrbayes_dir: Path, out_path: Path) -> None:
    """Генерировать HTML с филограммой на основе newick и клад."""
    # Парсим newick и строим простой rectangular phylogram
    try:
        from Bio import Phylo
        tree = Phylo.read(StringIO(newick), "newick")
    except Exception as e:
        print(f"  ⚠ Ошибка парсинга newick {group_key}: {e}", file=sys.stderr)
        return

    # Извлекаем листья и глубины
    tips = tree.get_terminals()
    n_tips = len(tips)
    if n_tips == 0:
        return

    # Сортируем листья в сюжетном порядке (DFS обход)
    leaf_order = []
    def dfs(clade):
        if clade.is_terminal():
            leaf_order.append(clade.name)
        else:
            for child in clade.clades:
                dfs(child)
    dfs(tree.root)

    # Высчитываем максимальную глубину (length от корня)
    def max_depth_of(clade, depth=0):
        if clade.is_terminal():
            return depth
        return max(max_depth_of(c, depth + (c.branch_length or 0)) for c in clade.clades) \
               if clade.clades else depth

    max_d = max_depth_of(tree.root)
    if max_d == 0:
        max_d = 1

    # Параметры layout
    left_margin = 40
    plot_w = 520
    scale = plot_w / max_d
    row_h = 56
    top = 60

    y_coords = {t: top + i*row_h for i, t in enumerate(leaf_order)}

    # Функция для вычисления X координаты узла
    def x_coord(clade, depth=0):
        return left_margin + depth*scale

    # Строим дерево с координатами
    def build_coords(clade, depth=0):
        x = x_coord(clade, depth)
        if clade.is_terminal():
            return {"name": clade.name, "x": x, "y": y_coords.get(clade.name, 0),
                    "type": "leaf", "confidence": clade.confidence}
        else:
            children = [build_coords(c, depth + (c.branch_length or 0)) for c in clade.clades]
            if children:
                avg_y = sum(c["y"] for c in children) / len(children)
            else:
                avg_y = 0
            return {"x": x, "y": avg_y, "type": "internal", "children": children,
                    "confidence": clade.confidence}

    tree_coords = build_coords(tree.root)

    # Определяем уверенные клады для подсветки
    confident_leaf_pairs = set()
    for clade_data in clades:
        if clade_data.get("posterior", 0) >= 0.95:
            leaves = tuple(sorted(clade_data["leaves"]))
            confident_leaf_pairs.add(leaves)

    # Генерируем SVG
    svg_lines = [
        f'<svg viewBox="0 0 760 {top + len(leaf_order)*row_h + 60}" role="img" aria-label="Дерево {group_key}">',
    ]

    # Рекурсивный рендер веток
    def render_branches(node, parent_x=None, parent_y=None, parent_confident=False):
        if node["type"] == "leaf":
            x, y = node["x"], node["y"]
            if parent_x is not None:
                # Горизонтальная ветка от родителя
                is_confident = False
                for pair in confident_leaf_pairs:
                    if node["name"] in pair:
                        is_confident = True
                        break
                css_class = "branch confident" if is_confident else "branch"
                svg_lines.append(f'<line class="{css_class}" x1="{parent_x}" y1="{parent_y}" '
                               f'x2="{x}" y2="{y}"></line>')
            # Точка листа
            css_class = "leaf-dot confident" if parent_confident else "leaf-dot"
            svg_lines.append(f'<circle class="{css_class}" cx="{x}" cy="{y}" r="3.5"></circle>')
            # Hover-подсказка
            svg_lines.append(f'<circle class="hit" cx="{x}" cy="{y}" r="10">'
                           f'<title>{node["name"]}</title></circle>')
            # Текстовая метка
            label_short = node["name"][:14] + ("…" if len(node["name"]) > 14 else "")
            svg_lines.append(f'<text class="leaf-label" x="{x+8}" y="{y-3}">{label_short}</text>')
        else:
            x, y = node["x"], node["y"]
            if parent_x is not None:
                # Вертикальная ветка от родителя к внутреннему узлу
                svg_lines.append(f'<line class="branch" x1="{parent_x}" y1="{parent_y}" '
                               f'x2="{x}" y2="{parent_y}"></line>')
                svg_lines.append(f'<line class="branch" x1="{x}" y1="{parent_y}" '
                               f'x2="{x}" y2="{y}"></line>')
            # Вертикальная линия между потомками
            if node["children"]:
                y_min = min(c["y"] for c in node["children"])
                y_max = max(c["y"] for c in node["children"])
                svg_lines.append(f'<line class="branch" x1="{x}" y1="{y_min}" '
                               f'x2="{x}" y2="{y_max}"></line>')
            # Узел с меткой
            if node["confidence"] is not None:
                pp = float(node["confidence"])
                if pp > 1:
                    pp = pp / 100
                pp = round(pp, 2)
                svg_lines.append(f'<text class="node-label" x="{x-5}" y="{y-5}" text-anchor="end">'
                               f'PP {pp}</text>')
            # Рекурсивный рендер потомков
            for child in node.get("children", []):
                child_is_conf = any(pair for pair in confident_leaf_pairs
                                  if set(pair) & set(get_leaves(child)))
                render_branches(child, x, y, child_is_conf)

    def get_leaves(node):
        if node["type"] == "leaf":
            return [node["name"]]
        return sum((get_leaves(c) for c in node.get("children", [])), [])

    render_branches(tree_coords)
    svg_lines.append('</svg>')

    # HTML шаблон
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
    background: var(--page);
    color: var(--text-primary);
    font-family: system-ui, -apple-system, "Segoe UI", sans-serif;
    padding: 28px 20px;
  }}
  .wrap {{ max-width: 860px; margin: 0 auto; }}
  .card {{
    background: var(--surface-1);
    border: 1px solid var(--border);
    border-radius: 14px;
    padding: 22px 24px;
    margin-bottom: 16px;
  }}
  h1 {{ font-size: 18px; margin: 0 0 4px; }}
  .sub {{
    color: var(--text-secondary);
    font-size: 13px;
    line-height: 1.5;
    margin-bottom: 4px;
  }}
  .meta {{
    display: flex; flex-wrap: wrap; gap: 8px 18px;
    font-size: 12.5px; color: var(--text-muted);
    margin-top: 10px; padding-top: 10px;
    border-top: 1px solid var(--gridline);
  }}
  .meta b {{ color: var(--text-secondary); font-weight: 600; }}
  .status-ok {{ color: var(--good); font-weight: 600; }}

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
        <span><b>уверенных клад</b> {len([c for c in clades if c.get("posterior", 0) >= 0.95])}</span>
      </div>
    </div>
    <div class="card">
      {"".join(svg_lines)}
    </div>
"""

    if clades:
        html += """    <div class="card">
      <h1 style="font-size:15px;">Уверенные клады (posterior ≥ 0.95)</h1>
      <table>
        <tr><th>листья</th><th class="n">posterior</th></tr>
"""
        for c in clades:
            if c.get("posterior", 0) >= 0.95:
                leaves_str = ", ".join(c["leaves"][:2]) + ("..." if len(c["leaves"]) > 2 else "")
                html += f'        <tr><td><span class="pill">✓</span> {leaves_str}</td>' \
                       f'<td class="n">{c.get("posterior", "—")}</td></tr>\n'
        html += "      </table>\n    </div>\n"

    html += "  </div>\n</div>"
    out_path.write_text(html, encoding="utf-8")


def run_test_cycle(n_groups: int = 3) -> None:
    """Запустить полный цикл: загрузка → mrbayes → confident_clades → визуализация."""
    print("=" * 70)
    print("ТЕСТОВЫЙ ЦИКЛ: nexus + MrBayes + уверенные клады")
    print("=" * 70)

    base_dir = Path(__file__).parent

    # 1. Загрузка тестовых групп
    print(f"\n1. Загрузка {n_groups} тестовых групп из BIOCAD.bigchallenges...")
    groups = download_test_groups(n_groups)
    if not groups:
        print("Ошибка: не удалось загрузить группы", file=sys.stderr)
        return

    # 2. Запуск mrbayes/run_mrbayes.py
    print(f"\n2. Запуск mrbayes/run_mrbayes.py (это займёт ~50s на группу)...")
    mrbayes_dir = base_dir / "mrbayes"
    result = subprocess.run(
        [sys.executable, str(base_dir / "mrbayes" / "run_mrbayes.py"),
         str(base_dir / "aligned_sequences"), "--out", str(mrbayes_dir),
         "--mb-ngen", "2000000"],
        capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(f"Ошибка mrbayes: {result.stderr}", file=sys.stderr)

    # 3. Запуск groups/confident_clades_report.py
    print(f"\n3. Запуск groups/confident_clades_report.py...")
    groups_dir = base_dir / "groups"
    result = subprocess.run(
        [sys.executable, str(base_dir / "groups" / "confident_clades_report.py"),
         "--mrbayes-dir", str(mrbayes_dir), "--posterior-min", "0.95",
         "--out", str(groups_dir / "report.json")],
        capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(f"Ошибка clades: {result.stderr}", file=sys.stderr)

    # 4. Генерация визуализаций
    print(f"\n4. Генерация HTML-визуализаций...")
    report = json.loads((groups_dir / "report.json").read_text(encoding="utf-8")) \
             if (groups_dir / "report.json").is_file() else {}

    for group_key in groups:
        con_tre = mrbayes_dir / f"{group_key}.nex.con.tre"
        names_tsv = mrbayes_dir / f"{group_key}.names.tsv"
        if not con_tre.is_file():
            print(f"  ⚠ {group_key}: .nex.con.tre не найден", file=sys.stderr)
            continue

        newick, _ = parse_nexus_tree(con_tre, names_tsv)
        if not newick:
            print(f"  ⚠ {group_key}: не удалось распарсить дерево", file=sys.stderr)
            continue

        clades = report.get(group_key, {}).get("mrbayes", {}).get("clades", [])
        out_html = mrbayes_dir / f"{group_key}.tree.html"
        generate_tree_html(group_key, newick, clades, mrbayes_dir, out_html)
        print(f"  ✓ {group_key}: {out_html}")

    # Итоги
    print("\n" + "=" * 70)
    print("РЕЗУЛЬТАТЫ:")
    print("=" * 70)
    for group_key in groups:
        con_tre = mrbayes_dir / f"{group_key}.nex.con.tre"
        report_file = groups_dir / "report.json"
        html_file = mrbayes_dir / f"{group_key}.tree.html"
        print(f"\n{group_key}:")
        if con_tre.is_file():
            print(f"  Дерево:       {con_tre.relative_to(Path.cwd())}")
        if html_file.is_file():
            print(f"  Визуализация: {html_file.relative_to(Path.cwd())}")
    if report_file.is_file():
        print(f"\nОтчёт всех групп: {report_file.relative_to(Path.cwd())}")
    print("\nОткройте .tree.html файлы в браузере для просмотра визуализаций.")


if __name__ == "__main__":
    import sys
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    run_test_cycle(n)
