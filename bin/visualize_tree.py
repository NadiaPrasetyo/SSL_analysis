#!/usr/bin/env python3
"""
Tree visualization tool for alifilter output
Generates a rectangular phylogenetic tree with coverage metrics and reference highlighting.
"""

import json
import os
import sys
import argparse
from pathlib import Path


def parse_newick(newick_str):
    """Parse Newick format tree string."""
    newick_str = newick_str.strip().rstrip(';')
    
    class Node:
        def __init__(self, name='', bootstrap=None):
            self.name = name
            self.bootstrap = bootstrap
            self.children = []
            self.distance = 0.0
    
    pos = 0
    
    def parse_node():
        nonlocal pos
        node = Node()
        
        if pos < len(newick_str) and newick_str[pos] == '(':
            pos += 1
            while pos < len(newick_str) and newick_str[pos] != ')':
                if newick_str[pos] == ',':
                    pos += 1
                else:
                    child = parse_node()
                    node.children.append(child)
            
            pos += 1
            
            label_start = pos
            while pos < len(newick_str) and newick_str[pos] not in '(),;:':
                pos += 1
            if pos > label_start:
                node.name = newick_str[label_start:pos]
        else:
            label_start = pos
            while pos < len(newick_str) and newick_str[pos] not in '(),;:':
                pos += 1
            node.name = newick_str[label_start:pos]
        
        if pos < len(newick_str) and newick_str[pos] == ':':
            pos += 1
            dist_start = pos
            while pos < len(newick_str) and newick_str[pos] not in '(),;':
                pos += 1
            try:
                node.distance = float(newick_str[dist_start:pos])
            except ValueError:
                node.distance = 0.0
        
        return node
    
    return parse_node()


def clean_seq_name(seq_name):
    """Remove leading sequence number (e.g., '861|ERR212785_...' -> 'ERR212785_...')"""
    if '|' in seq_name:
        parts = seq_name.split('|', 1)
        if parts[0].isdigit():
            return parts[1]
        return seq_name
    return seq_name


def extract_branch_from_seq(seq_name):
    """Extract branch (SSL3, SSL7, etc.) from sequence name"""
    clean_name = clean_seq_name(seq_name)
    if '|' in clean_name:
        return clean_name.split('|')[0]
    for prefix in ["SSL3", "SSL7", "SSL11", "SSL1", "SSL2"]:
        if clean_name.startswith(prefix):
            return prefix
    return "Unknown"


def collect_leaf_sequences(node):
    """Recursively collect all leaf sequences from a tree node"""
    if not node.children:
        return [clean_seq_name(node.name)] if node.name else []
    leaves = []
    for child in node.children:
        leaves.extend(collect_leaf_sequences(child))
    return leaves


def get_tree_depth(node):
    """Get maximum depth of tree from this node"""
    if not node.children:
        return node.distance
    return max(child.distance + get_tree_depth(child) for child in node.children)


def generate_rectangular_tree_svg(tree_path, coverage_json_path, svg_output_path):
    """Generate rectangular phylogenetic tree SVG."""
    
    # Hardcoded always-keep sequences
    always_keep_SSL3 = {
        "SSL3|CC1", "SSL3|CC5", "SSL3|CC8", "SSL3|CC22", "SSL3|CC30", "SSL3|CC93"
    }
    always_keep_SSL7 = {
        "SSL7|CC1", "SSL7|CC5", "SSL7|CC8", "SSL7|CC22", "SSL7|CC30", "SSL7|CC93"
    }
    always_keep_SSL11 = {
        "SSL11|CC1", "SSL11|CC5", "SSL11|CC8", "SSL11|CC22", "SSL11|CC30", "SSL11|CC93"
    }
    always_keep_set = always_keep_SSL3 | always_keep_SSL7 | always_keep_SSL11
    
    # Read tree and coverage data
    with open(tree_path, 'r') as f:
        newick_str = f.read().strip()
    tree = parse_newick(newick_str)
    
    with open(coverage_json_path, 'r') as f:
        sequence_coverage = json.load(f)
    
    # Build branch coverage cache
    branch_coverage_cache = {}
    for seq_name, seq_data in sequence_coverage.items():
        branch = extract_branch_from_seq(seq_name)
        if branch not in branch_coverage_cache:
            branch_coverage_cache[branch] = {'total': 0.0, 'count': 0}
        branch_coverage_cache[branch]['total'] += seq_data['coverage_percent']
        branch_coverage_cache[branch]['count'] += 1
    
    for branch in branch_coverage_cache:
        avg = (branch_coverage_cache[branch]['total'] / 
               branch_coverage_cache[branch]['count']) if branch_coverage_cache[branch]['count'] > 0 else 0.0
        branch_coverage_cache[branch] = avg
    
    # Color mapping
    branch_colors = {
        'SSL1': '#e74c3c',
        'SSL2': '#3498db',
        'SSL3': '#2980b9',
        'SSL5': '#27ae60',
        'SSL7': '#f39c12',
        'SSL11': '#d35400',
        'Unknown': '#95a5a6'
    }
    
    # Calculate tree layout
    all_leaves = collect_leaf_sequences(tree)
    num_leaves = len(all_leaves)
    
    # SVG dimensions
    margin = 100
    leaf_y_spacing = 40
    max_tree_distance = get_tree_depth(tree)
    x_scale = 400 / max(max_tree_distance, 0.1)
    
    svg_width = margin * 2 + max_tree_distance * x_scale
    svg_height = margin * 2 + num_leaves * leaf_y_spacing
    
    # Build rectangular tree layout
    leaf_counter = [0]
    node_positions = {}
    
    def layout_tree(node, depth=0):
        """Calculate positions for rectangular tree layout"""
        x = margin + depth * x_scale
        
        if not node.children:
            y = margin + leaf_counter[0] * leaf_y_spacing
            leaf_counter[0] += 1
            node_positions[id(node)] = (x, y, node.name)
            return x, y
        
        # Internal node: position between children
        child_positions = []
        for child in node.children:
            cx, cy = layout_tree(child, depth + child.distance)
            child_positions.append((cx, cy))
        
        # Y position is average of children
        y = sum(cy for cx, cy in child_positions) / len(child_positions)
        node_positions[id(node)] = (x, y, node.name)
        
        return x, y
    
    layout_tree(tree)
    
    # Generate SVG content
    svg_lines = [
        f'<svg version="1.1" xmlns="http://www.w3.org/2000/svg" width="{svg_width}" height="{svg_height}" viewBox="0 0 {svg_width} {svg_height}">',
        '<defs>',
        '<style>',
        '.branch { stroke: #000; stroke-width: 1.5; fill: none; }',
        '.leaf-circle { fill: #333; stroke: #000; stroke-width: 1; }',
        '.reference-highlight { fill: #FFE5B4; opacity: 0.6; }',
        '.leaf-label { font-family: Arial, sans-serif; font-size: 16px; fill: #000; }',
        '.reference-label { font-family: Arial, sans-serif; font-size: 16px; fill: #000; font-weight: bold; }',
        '.coverage-label { font-family: Arial, sans-serif; font-size: 12px; fill: #666; }',
        '</style>',
        '</defs>',
        '<rect width="100%" height="100%" fill="white"/>',
    ]
    
    # Draw tree edges
    def draw_tree(node):
        """Recursively draw tree branches"""
        if node not in node_positions:
            return
        
        x1, y1, _ = node_positions[id(node)]
        
        for child in node.children:
            if child not in node_positions:
                continue
            
            x2, y2, _ = node_positions[id(child)]
            
            # Get branch color from child
            clean_name = clean_seq_name(child.name)
            branch = extract_branch_from_seq(child.name)
            color = branch_colors.get(branch, '#95a5a6')
            
            # Rectangular connection: horizontal then vertical
            mid_x = x1 + (x2 - x1) * 0.7
            
            svg_lines.append(f'<line x1="{x1}" y1="{y1}" x2="{mid_x}" y2="{y1}" class="branch" stroke="{color}"/>')
            svg_lines.append(f'<line x1="{mid_x}" y1="{y1}" x2="{x2}" y2="{y2}" class="branch" stroke="{color}"/>')
            
            draw_tree(child)
    
    draw_tree(tree)
    
    # Draw leaf nodes
    for node_id, (x, y, name) in node_positions.items():
        if not isinstance(name, str) or not name:
            continue
        
        clean_name = clean_seq_name(name)
        branch = extract_branch_from_seq(name)
        color = branch_colors.get(branch, '#95a5a6')
        is_reference = clean_name in always_keep_set
        
        # Check if leaf (no children in original tree)
        for node in _get_all_nodes(tree):
            if id(node) == node_id and not node.children:
                # Draw highlight background for reference sequences
                if is_reference:
                    svg_lines.append(f'<rect x="{x-8}" y="{y-16}" width="24" height="24" class="reference-highlight" rx="4"/>')
                
                # Draw leaf circle
                svg_lines.append(f'<circle cx="{x}" cy="{y}" r="5" class="leaf-circle" fill="{color}"/>')
                
                # Draw label
                label = clean_name.split('|')[-1] if '|' in clean_name else clean_name
                label_class = 'reference-label' if is_reference else 'leaf-label'
                svg_lines.append(f'<text x="{x + 15}" y="{y + 6}" class="{label_class}">{label}</text>')
                
                # Add coverage info if available
                if clean_name in sequence_coverage:
                    cov = sequence_coverage[clean_name]['coverage_percent']
                    svg_lines.append(f'<text x="{x + 15}" y="{y + 22}" class="coverage-label">({cov:.1f}%)</text>')
                break
    
    svg_lines.append('</svg>')
    
    # Write SVG
    with open(svg_output_path, 'w') as f:
        f.write('\n'.join(svg_lines))
    
    return svg_output_path


def _get_all_nodes(node):
    """Helper to get all nodes in tree"""
    yield node
    for child in node.children:
        yield from _get_all_nodes(child)


def generate_html_visualization(coverage_json_path, svg_embed_path, html_output_path):
    """Generate HTML report with embedded SVG and coverage table."""
    
    # Hardcoded always-keep sequences
    always_keep_SSL3 = {
        "SSL3|CC1", "SSL3|CC5", "SSL3|CC8", "SSL3|CC22", "SSL3|CC30", "SSL3|CC93"
    }
    always_keep_SSL7 = {
        "SSL7|CC1", "SSL7|CC5", "SSL7|CC8", "SSL7|CC22", "SSL7|CC30", "SSL7|CC93"
    }
    always_keep_SSL11 = {
        "SSL11|CC1", "SSL11|CC5", "SSL11|CC8", "SSL11|CC22", "SSL11|CC30", "SSL11|CC93"
    }
    always_keep_set = always_keep_SSL3 | always_keep_SSL7 | always_keep_SSL11
    
    # Load coverage data
    with open(coverage_json_path, 'r') as f:
        sequence_coverage = json.load(f)
    
    # Create coverage table rows
    coverage_rows = ''
    for seq_name in sorted(sequence_coverage.keys()):
        stats = sequence_coverage[seq_name]
        num = stats['num']
        cov_pct = stats['coverage_percent']
        is_ref = seq_name in always_keep_set
        ref_marker = '⭐ ' if is_ref else ''
        
        coverage_rows += f'''
        <tr {'class="reference-row"' if is_ref else ''}>
            <td><strong>{ref_marker}{seq_name}</strong></td>
            <td class="numeric">{num}</td>
            <td class="numeric">{cov_pct:.2f}%</td>
        </tr>
        '''
    
    # Read SVG
    with open(svg_embed_path, 'r') as f:
        svg_content = f.read()
    
    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Phylogenetic Tree Visualization</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        @import url('https://fonts.googleapis.com/css2?family=Lato:wght@400;700&display=swap');
        
        body {{
            font-family: 'Lato', sans-serif;
            background: #ffffff;
            color: #333;
            padding: 20px;
            line-height: 1.6;
        }}
        
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
            padding: 30px;
        }}
        
        h1 {{
            text-align: center;
            margin-bottom: 8px;
            font-size: 2.2em;
            color: #2c3e50;
        }}
        
        .subtitle {{
            text-align: center;
            color: #7f8c8d;
            margin-bottom: 30px;
            font-size: 1em;
        }}
        
        .tree-section {{
            background: white;
            border: 2px solid #ecf0f1;
            border-radius: 8px;
            padding: 20px;
            margin-bottom: 30px;
            overflow-x: auto;
        }}
        
        .tree-section svg {{
            display: block;
            margin: 0 auto;
            background: white;
            width: 100%;
            height: auto;
        }}
        
        .coverage-section {{
            background: white;
            border: 2px solid #ecf0f1;
            border-radius: 8px;
            padding: 25px;
        }}
        
        .coverage-section h2 {{
            font-size: 1.5em;
            margin-bottom: 15px;
            color: #2c3e50;
        }}
        
        .coverage-section > p {{
            margin-bottom: 20px;
            color: #555;
            font-size: 0.95em;
        }}
        
        .coverage-table {{
            width: 100%;
            border-collapse: collapse;
            margin-bottom: 20px;
        }}
        
        .coverage-table thead {{
            background: #ecf0f1;
            border-bottom: 2px solid #bdc3c7;
        }}
        
        .coverage-table th {{
            padding: 12px;
            text-align: left;
            font-weight: 700;
            color: #2c3e50;
            font-size: 0.9em;
        }}
        
        .coverage-table td {{
            padding: 11px 12px;
            border-bottom: 1px solid #ecf0f1;
            font-size: 0.95em;
        }}
        
        .coverage-table td.numeric {{
            text-align: right;
            font-family: 'Courier New', monospace;
        }}
        
        .coverage-table tbody tr:hover {{
            background: #f5f5f5;
        }}
        
        .coverage-table tbody tr.reference-row {{
            background: #FFF8DC;
        }}
        
        .coverage-table tbody tr.reference-row:hover {{
            background: #FFE5B4;
        }}
        
        .coverage-explanation {{
            background: #f9f9f9;
            border-left: 4px solid #3498db;
            padding: 15px;
            margin-top: 20px;
            border-radius: 4px;
            font-size: 0.9em;
            color: #555;
        }}
        
        .coverage-explanation h3 {{
            color: #2c3e50;
            margin-bottom: 10px;
            font-size: 1em;
        }}
        
        .explanation-item {{
            margin: 8px 0;
            line-height: 1.5;
        }}
        
        .legend {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 15px;
            margin-top: 25px;
            padding-top: 20px;
            border-top: 1px solid #ecf0f1;
        }}
        
        .legend-item {{
            display: flex;
            align-items: center;
            gap: 10px;
            font-size: 0.95em;
        }}
        
        .legend-color {{
            width: 18px;
            height: 18px;
            border-radius: 50%;
            border: 1px solid #333;
        }}
        
        .legend-item.highlight {{
            display: flex;
            align-items: center;
            gap: 10px;
        }}
        
        .reference-box {{
            width: 18px;
            height: 18px;
            background: #FFE5B4;
            border-radius: 3px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 Phylogenetic Tree Analysis</h1>
        <div class="subtitle">Rectangular tree layout with sequence coverage metrics</div>
        
        <div class="tree-section">
            {svg_content}
        </div>
        
        <div class="coverage-section">
            <h2>Sequence Coverage Report</h2>
            <p>Coverage metric: Number of sequences this representative accounts for, divided by total sequences in the original alignment.</p>
            
            <table class="coverage-table">
                <thead>
                    <tr>
                        <th>Sequence Name</th>
                        <th>Represents</th>
                        <th>Coverage %</th>
                    </tr>
                </thead>
                <tbody>
                    {coverage_rows}
                </tbody>
            </table>
            
            <div class="coverage-explanation">
                <h3>📊 Understanding the Metrics</h3>
                <div class="explanation-item">
                    <strong>Represents:</strong> The count of sequences this representative accounts for (1 = itself, plus any sequences merged with it due to >90% identity)
                </div>
                <div class="explanation-item">
                    <strong>Coverage %:</strong> (Sequences Represented / Total Original Sequences) × 100
                </div>
                <div class="explanation-item" style="margin-top: 12px; padding-top: 12px; border-top: 1px solid #e0e0e0;">
                    <strong>Example:</strong> If a representative accounts for 25 sequences out of 1000 total, its coverage is 2.5%
                </div>
                <div class="explanation-item" style="color: #d4af37; margin-top: 8px;">
                    <strong>⭐ = Always-Keep Reference Sequence</strong> — Retained regardless of redundancy
                </div>
            </div>
            
            <div class="legend">
                <div class="legend-item">
                    <div class="legend-color" style="background: #2980b9;"></div>
                    <span>SSL3</span>
                </div>
                <div class="legend-item">
                    <div class="legend-color" style="background: #f39c12;"></div>
                    <span>SSL7</span>
                </div>
                <div class="legend-item">
                    <div class="legend-color" style="background: #d35400;"></div>
                    <span>SSL11</span>
                </div>
                <div class="legend-item">
                    <div class="legend-color" style="background: #27ae60;"></div>
                    <span>SSL5</span>
                </div>
                <div class="legend-item highlight">
                    <div class="reference-box"></div>
                    <span>Reference Highlight</span>
                </div>
            </div>
        </div>
    </div>
</body>
</html>
'''
    
    with open(html_output_path, 'w') as f:
        f.write(html)
    
    print(f"✅ HTML report saved to: {html_output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate rectangular tree visualization from alifilter output'
    )
    parser.add_argument('--tree', required=True, help='Path to Newick tree file')
    parser.add_argument('--coverage', required=True, help='Path to sequence coverage JSON file')
    parser.add_argument('--svg-output', required=True, help='Output SVG file path')
    parser.add_argument('--html-output', required=True, help='Output HTML file path')
    
    args = parser.parse_args()
    
    # Generate SVG
    svg_path = generate_rectangular_tree_svg(args.tree, args.coverage, args.svg_output)
    print(f"✅ SVG saved to: {svg_path}")
    
    # Generate HTML
    generate_html_visualization(args.coverage, svg_path, args.html_output)


if __name__ == '__main__':
    main()