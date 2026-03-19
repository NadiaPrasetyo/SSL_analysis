"""
Heatmap with phylogenetic tree (Newick) plotted on the left.
- Y-axis: accessions (ordered by tree leaf order)
- X-axis: features
- Colour: z-score per feature (blue=-1, white=0, red=+1)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm
from scipy.cluster.hierarchy import dendrogram
from scipy.spatial.distance import squareform
import argparse

def main(input_csv, newick_tree_file, output_file):

    #  DATA
    df = pd.read_csv(input_csv)
    df = df.sort_values(by=["accession", "feature"], ignore_index=True)

    # ── 2.  NEWICK TREE  ──────────────────────────────────────────────────────────
    with open(newick_tree_file) as f:
        newick = f.read()

    # ── 3.  MINIMAL NEWICK PARSER  ────────────────────────────────────────────────
    class Node:
        def __init__(self, name="", length=0.0):
            self.name = name
            self.length = length
            self.children = []
            self.x = 0.0
            self.y = 0.0

    def parse_newick(s):
        s = s.strip().rstrip(";")
        def helper(s):
            node = Node()
            if s.startswith("("):
                depth, i = 0, 0
                for i, c in enumerate(s):
                    if c == "(": depth += 1
                    elif c == ")": depth -= 1
                    if depth == 0: break
                inner = s[1:i]
                rest  = s[i+1:]
                parts, depth2, start = [], 0, 0
                for j, c in enumerate(inner):
                    if c == "(": depth2 += 1
                    elif c == ")": depth2 -= 1
                    elif c == "," and depth2 == 0:
                        parts.append(inner[start:j])
                        start = j + 1
                parts.append(inner[start:])
                node.children = [helper(p) for p in parts]
                if ":" in rest:
                    label, length = rest.split(":", 1)
                    node.name = label.strip()
                    node.length = float(length.strip())
                else:
                    node.name = rest.strip()
            else:
                if ":" in s:
                    name, length = s.split(":", 1)
                    node.name = name.strip()
                    node.length = float(length.strip())
                else:
                    node.name = s.strip()
            return node
        return helper(s)

    root = parse_newick(newick)

    def assign_x(node, parent_x=0.0):
        node.x = parent_x + node.length
        for c in node.children:
            assign_x(c, node.x)

    assign_x(root)

    leaf_counter = [0]
    leaf_order = []

    def assign_y(node):
        if not node.children:
            node.y = leaf_counter[0]
            leaf_order.append(node.name)
            leaf_counter[0] += 1
        else:
            for c in node.children:
                assign_y(c)
            node.y = np.mean([c.y for c in node.children])

    assign_y(root)

    n_leaves = leaf_counter[0]

    # ── 4.  PIVOT & Z-SCORE  ─────────────────────────────────────────────────────
    pivot = df.pivot_table(index="accession", columns="feature", values="score")
    pivot = pivot.loc[pivot.index.isin(leaf_order)]   # drop accessions not in tree
    pivot = pivot.loc[leaf_order]                      # reorder rows to match tree

    from scipy.stats import zscore as _zscore
    z = pivot.apply(lambda col: _zscore(col, nan_policy="omit"), axis=0)
    z = z.clip(-1, 1)   # clamp to [-1, 1]

    features = list(z.columns)
    n_acc    = len(leaf_order)
    n_feat   = len(features)

    # ── 5.  FIGURE LAYOUT  ───────────────────────────────────────────────────────
    cell_h   = 0.45
    cell_w   = 0.70
    tree_w   = 3.0
    cbar_w   = 0.4
    pad      = 0.3

    fig_h = max(6, n_acc * cell_h + 2.0)
    fig_w = tree_w + n_feat * cell_w + cbar_w + pad * 3

    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")

    gs = fig.add_gridspec(
        1, 3,
        width_ratios=[tree_w, n_feat * cell_w, cbar_w],
        left=0.08, right=0.97,
        top=0.95, bottom=0.22,
        wspace=0.03,
    )

    ax_tree = fig.add_subplot(gs[0])
    ax_heat = fig.add_subplot(gs[1])
    ax_cbar = fig.add_subplot(gs[2])

    # ── 6.  DRAW TREE  ────────────────────────────────────────────────────────────
    tree_color = "black"
    branch_lw  = 1.3

    def draw_tree_proper(ax, root):
        def _draw(node):
            for c in node.children:
                ax.plot([node.x, c.x], [c.y, c.y],
                        color=tree_color, lw=branch_lw, solid_capstyle="round")
                _draw(c)
            if node.children:
                y_vals = [c.y for c in node.children]
                ax.plot([node.x, node.x], [min(y_vals), max(y_vals)],
                        color=tree_color, lw=branch_lw, solid_capstyle="round")
        _draw(root)

    draw_tree_proper(ax_tree, root)

    def all_nodes(node):
        yield node
        for c in node.children:
            yield from all_nodes(c)

    max_x = max(n.x for n in all_nodes(root))

    ax_tree.set_xlim(-max_x * 0.05, max_x * 1.6)   # extra right space for labels
    ax_tree.set_ylim(-0.5, n_leaves - 0.5)
    ax_tree.set_facecolor("white")
    ax_tree.axis("off")

    # leaf labels on the right of the tree (between tree and heatmap)
    group_colors = {"SSL3": "#d62728", "SSL7": "#2ca02c", "SSL11": "#1f77b4"}

    for leaf_name in leaf_order:
        idx = leaf_order.index(leaf_name)
        grp = leaf_name.split("|")[0]
        col = group_colors.get(grp, "black")
        ax_tree.text(max_x * 1.03, idx, leaf_name,
                    va="center", ha="left", fontsize=7.5,
                    color=col, fontfamily="monospace")

    # ── 7.  DRAW HEATMAP  ─────────────────────────────────────────────────────────
    norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    cmap = plt.cm.bwr

    z_arr = z.values

    im = ax_heat.imshow(
        z_arr,
        aspect="auto",
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
        origin="upper",
    )

    for i in range(n_acc + 1):
        ax_heat.axhline(i - 0.5, color="white", lw=0.5)
    for j in range(n_feat + 1):
        ax_heat.axvline(j - 0.5, color="white", lw=0.5)

    ax_heat.set_facecolor("white")
    ax_heat.set_xticks(range(n_feat))
    ax_heat.set_xticklabels(features, rotation=45, ha="right", fontsize=7.5,
                            color="black")
    ax_heat.xaxis.set_label_position("bottom")
    ax_heat.xaxis.tick_bottom()
    ax_heat.set_yticks([])
    ax_heat.tick_params(axis="x", length=0, pad=3)

    # ── 8.  COLORBAR  ────────────────────────────────────────────────────────────
    cb = fig.colorbar(im, cax=ax_cbar)
    cb.set_label("Z-score", color="black", fontsize=8, labelpad=6)
    cb.ax.yaxis.set_tick_params(color="black", labelsize=7)
    plt.setp(plt.getp(cb.ax.axes, "yticklabels"), color="black")
    cb.ax.set_facecolor("white")
    cb.outline.set_edgecolor("#aaa")

    # ── 9.  LEGEND & TITLE  ──────────────────────────────────────────────────────
    legend_handles = [
        mpatches.Patch(color=v, label=k) for k, v in group_colors.items()
    ]
    fig.legend(handles=legend_handles, loc="lower left", framealpha=0.8,
            labelcolor="black", fontsize=8,
            bbox_to_anchor=(0.01, 0.01), ncol=3,
            facecolor="white", edgecolor="#aaa")

    fig.suptitle("Accession Feature Heatmap  (Z-score per feature)",
                color="black", fontsize=13, fontweight="bold", y=1)

    fig.patch.set_facecolor("white")

    plt.savefig(output_file,
                dpi=180, bbox_inches="tight", facecolor=fig.get_facecolor())
    print(f"Saved → {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_csv", type=str, required=True, help="input CSV file")
    parser.add_argument("--tree", type=str, required=True, help="input tree file (Newick format)")
    parser.add_argument("-o", "--output_file", type=str, default="heatmap_tree.png", help="output file name (default: heatmap_tree.png)")
    args = parser.parse_args()

    main(args.input_csv, args.tree, args.output_file)