#!/usr/bin/env python3
"""
Tree visualization tool for alifilter output
Generates a rectangular phylogenetic tree with proper Newick parsing and layout.
"""

import json
import os
import argparse
import matplotlib.pyplot as plt


# ---------------------------
# Newick parsing
# ---------------------------

class Node:
    def __init__(self, name="", length=0.0):
        self.name = name
        self.length = length
        self.children = []
        self.parent = None


def parse_newick(newick):
    """Parse Newick string into a tree of Node objects."""

    tokens = []
    buf = ""
    for c in newick.strip():
        if c in "(),:;":
            if buf:
                tokens.append(buf)
                buf = ""
            tokens.append(c)
        else:
            buf += c
    if buf:
        tokens.append(buf)

    stack = []
    current = Node()

    i = 0
    while i < len(tokens):
        t = tokens[i]

        if t == "(":
            child = Node()
            child.parent = current
            current.children.append(child)
            stack.append(current)
            current = child

        elif t == ",":
            parent = stack[-1]
            sibling = Node()
            sibling.parent = parent
            parent.children.append(sibling)
            current = sibling

        elif t == ")":
            current = stack.pop()

        elif t == ":":
            i += 1
            current.length = float(tokens[i])

        elif t == ";":
            pass

        else:
            current.name = t

        i += 1

    return current


# ---------------------------
# Helpers
# ---------------------------

def clean_seq_name(seq_name):
    if "|" in seq_name:
        parts = seq_name.split("|", 1)
        if parts[0].isdigit():
            return parts[1]
    return seq_name


def extract_branch_from_seq(seq_name):
    name = clean_seq_name(seq_name)
    if "|" in name:
        return name.split("|")[0]
    return "Unknown"


def assign_y_positions(node, y=0, step=1, positions=None):
    """Assign vertical positions to leaves."""
    if positions is None:
        positions = {}

    if not node.children:
        positions[node] = y
        return y + step, positions

    for child in node.children:
        y, positions = assign_y_positions(child, y, step, positions)

    positions[node] = sum(positions[c] for c in node.children) / len(node.children)
    return y, positions


def assign_x_positions(node, parent_x=0, positions=None):
    """Assign horizontal positions using branch lengths."""
    if positions is None:
        positions = {}

    x = parent_x + node.length
    positions[node] = x

    for child in node.children:
        assign_x_positions(child, x, positions)

    return positions


# ---------------------------
# Main plotting
# ---------------------------

def generate_rectangular_tree_svg(tree_path, coverage_json_path, svg_output_path):

    with open(tree_path) as f:
        newick = f.read().strip()

    root = parse_newick(newick)

    with open(coverage_json_path) as f:
        coverage_data = json.load(f)

    # Normalize coverage lookup
    coverage_lookup = {
        clean_seq_name(k): v["coverage_percent"]
        for k, v in coverage_data.items()
    }


    highlight_colors = {
        'SSL3': '#4A90E2',
        'SSL7': '#7ED321',
        'SSL11': '#F5D547',
        'Unknown': '#000000'
    }

    # Layout
    _, y_pos = assign_y_positions(root)
    x_pos = assign_x_positions(root)

    fig, ax = plt.subplots(figsize=(10, 12))

    def draw(node):
        for child in node.children:
            # horizontal line
            ax.plot(
                [x_pos[node], x_pos[child]],
                [y_pos[child], y_pos[child]],
                color="black"
            )

            # vertical connector
            ax.plot(
                [x_pos[node], x_pos[node]],
                [y_pos[node], y_pos[child]],
                color="black"
            )

            draw(child)

    draw(root)

    # Labels
    def draw_labels(node):
        if not node.children:
            clean_name = clean_seq_name(node.name)
            branch = extract_branch_from_seq(node.name)

            label = clean_name

            if clean_name in coverage_lookup:
                cov = coverage_lookup[clean_name]
                label += f" ({cov:.2f}%)"

            color = highlight_colors.get(branch, '#000000')

            ax.text(
                x_pos[node] + 0.01,
                y_pos[node],
                label,
                verticalalignment='center',
                fontsize=8,
                color=color
            )
        else:
            for c in node.children:
                draw_labels(c)

    draw_labels(root)

    ax.set_ylim(-2, max(y_pos.values()) + 1)
    ax.set_xlim(0, max(x_pos.values()) * 1.1)

    ax.axis("off")

    # --- Tree scale bar ---
    max_x = max(x_pos.values())
    # Choose a round scale length: 1 if tree is small, 10 if large
    scale_length = 10 if max_x > 50 else 1

    # Position in the bottom-left
    bar_x = max_x * 0.02
    bar_y = -1.5  # just inside the bottom y-limit

    ax.plot(
        [bar_x, bar_x + scale_length],
        [bar_y, bar_y],
        color="black",
        linewidth=1.5,
        solid_capstyle="butt"
    )
    # Tick marks at each end
    tick_height = 0.2
    ax.plot([bar_x, bar_x], [bar_y - tick_height, bar_y + tick_height],
            color="black", linewidth=1.5)
    ax.plot([bar_x + scale_length, bar_x + scale_length],
            [bar_y - tick_height, bar_y + tick_height],
            color="black", linewidth=1.5)
    # Label centred above the bar
    # Label with legend text to the left of the bar
    ax.text(
        bar_x,
        bar_y + 0.35,
        f"Tree scale: {scale_length}",
        ha="left",
        va="bottom",
        fontsize=8
    )

    os.makedirs(os.path.dirname(svg_output_path), exist_ok=True)
    plt.savefig(svg_output_path, format="svg", bbox_inches="tight")
    plt.close()

    print(f"Saved SVG to {svg_output_path}")


# ---------------------------
# CLI
# ---------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Generate rectangular tree visualization from alifilter output'
    )
    parser.add_argument('--tree', required=True)
    parser.add_argument('--coverage', required=True)
    parser.add_argument('-o', '--output', required=True)

    args = parser.parse_args()

    generate_rectangular_tree_svg(
        args.tree,
        args.coverage,
        args.output
    )


if __name__ == '__main__':
    main()