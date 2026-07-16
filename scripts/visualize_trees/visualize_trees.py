import io
import os
import re
import sys

import toyplot.html
import toytree

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>
  body {{ margin: 0; background: white; }}
  .toyplot svg {{ background: white !important; }}
</style>
</head>
<body>
{content}
</body>
</html>"""

if len(sys.argv) != 3:
    print("Usage: python visualize_trees.py <input_dir> <output_dir>")
    sys.exit(1)

INPUT_DIR = sys.argv[1]
OUTPUT_DIR = sys.argv[2]

if not os.path.isdir(INPUT_DIR):
    print(f"Error: input directory {INPUT_DIR} not found")
    sys.exit(1)


def parse_iqtree_support(name):
    """Parse IQ-TREE support values (UFBoot/SH-aLRT format: '95/80')."""
    if not name or "/" not in str(name):
        return None, None
    parts = str(name).split("/")
    try:
        ufboot = float(parts[0])
        sh_alrt = float(parts[1])
        return ufboot, sh_alrt
    except (ValueError, IndexError):
        return None, None


def parse_mrbayes_support(name):
    """Parse MrBayes posterior probability from node label.
    MrBayes formats: '[posterior=0.95]' or just '0.95'"""
    if not name:
        return None
    name_str = str(name)
    # Match [posterior=0.95] or [0.95] or just 0.95
    match = re.search(r'\[?posterior\s*=\s*([0-9.]+)\]?|\[?([0-9.]+)\]?', name_str)
    if match:
        try:
            return float(match.group(1) or match.group(2))
        except (ValueError, IndexError):
            pass
    return None


def get_tree_type(filename):
    """Determine tree type from filename."""
    fname = filename.lower()
    if fname.endswith(".treefile"):
        return "iqtree"
    # MrBayes common extensions
    if fname.endswith((".con.tre", ".t", ".nex", ".tre", ".tree")):
        return "mrbayes"
    return None


def parse_support_for_tree(tree, tree_type):
    """Parse node support values based on tree type."""
    names = tree.get_node_data("name").to_list()
    if tree_type == "iqtree":
        return [parse_iqtree_support(n) for n in names]
    elif tree_type == "mrbayes":
        return [(parse_mrbayes_support(n), None) for n in names]
    return [(None, None) for _ in names]


found = 0
for root, dirs, files in os.walk(INPUT_DIR):
    for f in files:
        tree_type = get_tree_type(f)
        if tree_type is None:
            continue
        found = 1
        tree_path = os.path.join(root, f)
        rel = os.path.relpath(root, INPUT_DIR)
        outdir = os.path.join(OUTPUT_DIR, rel)
        os.makedirs(outdir, exist_ok=True)
        basename = f.replace(".treefile", "").replace(".con.tre", "").replace(".tre", "").replace(".t", "").replace(".nex", "")

        tree = toytree.tree(tree_path).ladderize()

        nleaves = tree.ntips
        height = max(400, nleaves * 35)

        names = tree.get_node_data("name").to_list()
        dists = tree.get_node_data("dist").to_list()

        node_labels = []
        node_sizes = []
        node_colors = []

        supports = parse_support_for_tree(tree, tree_type)

        for n, d, support in zip(names, dists, supports):
            if tree_type == "iqtree":
                ufboot, sh_alrt = support
                if ufboot is not None:
                    avg_support = (ufboot + sh_alrt) / 2
                    node_labels.append(f"UF: {ufboot:.0f} / SH: {sh_alrt:.0f}")
                    node_sizes.append(10 + avg_support * 0.2)
                    if avg_support >= 80:
                        node_colors.append("#1b7837")
                    elif avg_support >= 50:
                        node_colors.append("#a6d96a")
                    else:
                        node_colors.append("#ca0020")
                else:
                    node_labels.append("")
                    node_sizes.append(5)
                    node_colors.append("white")
            elif tree_type == "mrbayes":
                posterior, _ = support
                if posterior is not None:
                    pp_pct = posterior * 100
                    node_labels.append(f"PP: {posterior:.2f}")
                    node_sizes.append(10 + pp_pct * 0.2)
                    if pp_pct >= 95:
                        node_colors.append("#1b7837")
                    elif pp_pct >= 75:
                        node_colors.append("#a6d96a")
                    else:
                        node_colors.append("#ca0020")
                else:
                    node_labels.append("")
                    node_sizes.append(5)
                    node_colors.append("white")

        canvas, axes, mark = tree.draw(
            width=1200,
            height=height,
            tip_labels_align=True,
            tip_labels_style={"font-size": "11px", "fill": "#333"},
            node_labels=node_labels,
            node_sizes=node_sizes,
            node_colors=node_colors,
            node_labels_style={"font-size": "8px", "fill": "#222"},
            edge_style={"stroke": "#555", "stroke-width": 2},
        )

        canvas.background = "white"

        buf = io.BytesIO()
        toyplot.html.render(canvas, buf)
        content = buf.getvalue().decode("utf-8")

        html_path = os.path.join(outdir, f"{basename}.html")
        with open(html_path, "w") as out:
            out.write(HTML_TEMPLATE.format(title=basename, content=content))
        print(f"  Saved {html_path}")

if not found:
    print(f"No supported tree files found in {INPUT_DIR} (.treefile, .con.tre, .t, .nex, .tre)")
    sys.exit(1)
