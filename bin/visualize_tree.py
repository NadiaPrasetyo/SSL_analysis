#!/usr/bin/env python3
"""
Tree visualization tool for alifilter output
Generates a rectangular phylogenetic tree with proper spacing and reference highlighting.
"""

import json
import os
import argparse


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
    """Remove leading sequence number."""
    if '|' in seq_name:
        parts = seq_name.split('|', 1)
        if parts[0].isdigit():
            return parts[1]
        return seq_name
    return seq_name


def extract_branch_from_seq(seq_name):
    """Extract branch type from sequence name."""
    clean_name = clean_seq_name(seq_name)
    if '|' in clean_name:
        return clean_name.split('|')[0]
    for prefix in ["SSL3", "SSL7", "SSL11", "SSL1", "SSL2"]:
        if clean_name.startswith(prefix):
            return prefix
    return "Unknown"


def collect_leaf_count(node):
    """Count number of leaves below this node."""
    if not node.children:
        return 1
    return sum(collect_leaf_count(child) for child in node.children)


def get_max_distance(node):
    """Get maximum branch distance in tree."""
    if not node.children:
        return 0
    return max(node.distance + get_max_distance(child) for child in node.children)


def generate_rectangular_tree_svg(tree_path, coverage_json_path, svg_output_path):
    """Generate rectangular phylogenetic tree SVG with proper spacing."""
    
    # Always-keep sequences
    always_keep_set = {
        "SSL3|CC1", "SSL3|CC5", "SSL3|CC8", "SSL3|CC22", "SSL3|CC30", "SSL3|CC93",
        "SSL7|CC1", "SSL7|CC5", "SSL7|CC8", "SSL7|CC22", "SSL7|CC30", "SSL7|CC93",
        "SSL11|CC1", "SSL11|CC5", "SSL11|CC8", "SSL11|CC22", "SSL11|CC30", "SSL11|CC93"
    }
    
    # Read input files
    with open(tree_path, 'r') as f:
        tree = parse_newick(f.read().strip())
    
    with open(coverage_json_path, 'r') as f:
        sequence_coverage = json.load(f)
    
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
    
    highlight_colors = {
        'SSL3': '#4A90E2',
        'SSL7': '#7ED321',
        'SSL11': '#F5D547',
        'SSL1': '#FF6B6B',
        'SSL2': '#4ECDC4',
        'SSL5': '#95E1D3',
        'Unknown': '#CCCCCC'
    }
    
    # Layout parameters - make SVG VERY WIDE to show branch distances
    margin_left = 80
    margin_top = 40
    leaf_height = 30
    x_scale = 50  # Much larger for branch distance visibility
    label_right_margin = 300
    
    # Calculate dimensions
    total_leaves = collect_leaf_count(tree)
    max_distance = get_max_distance(tree)
    
    svg_height = margin_top * 2 + total_leaves * leaf_height
    svg_width = margin_left + (max_distance * x_scale) + label_right_margin
    
    # Layout tree
    leaf_y = [margin_top]
    node_positions = {}
    
    def layout_tree(node, depth=0):
        """Calculate node positions."""
        x = margin_left + depth * x_scale
        
        if not node.children:
            y = leaf_y[0]
            leaf_y[0] += leaf_height
            node_positions[id(node)] = (x, y)
            return x, y
        
        # Internal node at midpoint of children
        child_positions = []
        for child in node.children:
            cx, cy = layout_tree(child, depth + node.distance)
            child_positions.append((cx, cy))
        
        avg_y = sum(cy for cx, cy in child_positions) / len(child_positions)
        x_internal = margin_left + (depth + node.distance) * x_scale
        node_positions[id(node)] = (x_internal, avg_y)
        
        return x_internal, avg_y
    
    layout_tree(tree)
    
    # Build SVG
    svg = [
        f'<svg version="1.1" xmlns="http://www.w3.org/2000/svg" width="{svg_width:.0f}" height="{svg_height:.0f}" viewBox="0 0 {svg_width:.0f} {svg_height:.0f}">',
        '<defs><style>',
        '.branch { stroke: #333; stroke-width: 2; fill: none; }',
        '.leaf-circle { stroke: #000; stroke-width: 1.5; }',
        '.leaf-label { font-family: Arial, sans-serif; font-size: 16px; fill: #000; }',
        '.reference-label { font-family: Arial, sans-serif; font-size: 16px; fill: #fff; font-weight: bold; }',
        '.coverage-label { font-family: Arial, sans-serif; font-size: 12px; fill: #666; }',
        '</style></defs>',
        '<rect width="100%" height="100%" fill="white"/>',
    ]
    
    # Draw branches
    def draw_branches(node):
        """Draw rectangular tree branches."""
        if id(node) not in node_positions:
            return
        
        x1, y1 = node_positions[id(node)]
        
        for child in node.children:
            if id(child) not in node_positions:
                continue
            
            x2, y2 = node_positions[id(child)]
            
            # Horizontal line from parent to child x-position
            svg.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y1:.1f}" class="branch"/>')
            # Vertical line from parent y to child y
            svg.append(f'<line x1="{x2:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" class="branch"/>')
            
            draw_branches(child)
    
    draw_branches(tree)
    
    # Draw leaf nodes
    def draw_leaves(node):
        """Draw leaf nodes and labels with highlighting."""
        if node.children:
            for child in node.children:
                draw_leaves(child)
            return
        
        if id(node) not in node_positions:
            return
        
        x, y = node_positions[id(node)]
        clean_name = clean_seq_name(node.name)
        branch = extract_branch_from_seq(node.name)
        is_reference = clean_name in always_keep_set
        
        # Draw circle
        color = branch_colors.get(branch, '#95a5a6')
        svg.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5" class="leaf-circle" fill="{color}"/>')
        
        # Draw label (clean version)
        label = clean_name.split('|')[-1] if '|' in clean_name else clean_name
        
        if is_reference:
            # Highlighted background for reference sequences
            highlight_color = highlight_colors.get(branch, '#CCCCCC')
            label_width = len(label) * 9  # Approximate width
            svg.append(f'<rect x="{x + 12}" y="{y - 18}" width="{label_width + 8}" height="22" fill="{highlight_color}" rx="3"/>')
            svg.append(f'<text x="{x + 16}" y="{y + 3}" class="reference-label">{label}</text>')
        else:
            svg.append(f'<text x="{x + 12}" y="{y + 5}" class="leaf-label">{label}</text>')
        
        # Add coverage info
        if clean_name in sequence_coverage:
            cov = sequence_coverage[clean_name]['coverage_percent']
            svg.append(f'<text x="{x + 12}" y="{y + 20}" class="coverage-label">({cov:.1f}%)</text>')
    
    draw_leaves(tree)
    
    svg.append('</svg>')
    
    # Write file
    with open(svg_output_path, 'w') as f:
        f.write('\n'.join(svg))
    
    print(f"✅ SVG saved to: {svg_output_path}")


def generate_html_visualization(coverage_json_path, svg_output_path, html_output_path):
    """Generate HTML report."""
    
    always_keep_set = {
        "SSL3|CC1", "SSL3|CC5", "SSL3|CC8", "SSL3|CC22", "SSL3|CC30", "SSL3|CC93",
        "SSL7|CC1", "SSL7|CC5", "SSL7|CC8", "SSL7|CC22", "SSL7|CC30", "SSL7|CC93",
        "SSL11|CC1", "SSL11|CC5", "SSL11|CC8", "SSL11|CC22", "SSL11|CC30", "SSL11|CC93"
    }
    
    with open(coverage_json_path, 'r') as f:
        sequence_coverage = json.load(f)
    
    # Create table rows
    coverage_rows = ''
    for seq_name in sorted(sequence_coverage.keys()):
        stats = sequence_coverage[seq_name]
        num = stats['num']
        cov_pct = stats['coverage_percent']
        is_ref = seq_name in always_keep_set
        ref_marker = '⭐' if is_ref else ''
        
        coverage_rows += f'''
        <tr {'class="reference-row"' if is_ref else ''}>
            <td><strong>{ref_marker} {seq_name}</strong></td>
            <td class="numeric">{num}</td>
            <td class="numeric">{cov_pct:.2f}%</td>
        </tr>
        '''
    
    # Read SVG
    with open(svg_output_path, 'r') as f:
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
            max-width: 100%;
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
            max-height: 1200px;
            overflow-y: auto;
        }}
        
        .tree-section svg {{
            display: block;
            background: white;
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
            font-size: 0.95em;
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
        }}
        
        .coverage-table td {{
            padding: 11px 12px;
            border-bottom: 1px solid #ecf0f1;
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
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin-top: 25px;
            padding-top: 20px;
            border-top: 1px solid #ecf0f1;
            font-size: 0.95em;
        }}
        
        .legend-item {{
            display: flex;
            align-items: center;
            gap: 10px;
        }}
        
        .legend-color {{
            width: 18px;
            height: 18px;
            border-radius: 50%;
            border: 1px solid #333;
            flex-shrink: 0;
        }}
        
        .reference-box {{
            width: 18px;
            height: 18px;
            border-radius: 2px;
            border: 1px solid #333;
            flex-shrink: 0;
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
                <div class="legend-item">
                    <div class="reference-box" style="background: #4A90E2; border: none;"></div>
                    <span>SSL3 Reference (Highlighted)</span>
                </div>
                <div class="legend-item">
                    <div class="reference-box" style="background: #7ED321; border: none;"></div>
                    <span>SSL7 Reference (Highlighted)</span>
                </div>
                <div class="legend-item">
                    <div class="reference-box" style="background: #F5D547; border: none;"></div>
                    <span>SSL11 Reference (Highlighted)</span>
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
    generate_rectangular_tree_svg(args.tree, args.coverage, args.svg_output)
    
    # Generate HTML
    generate_html_visualization(args.coverage, args.svg_output, args.html_output)


if __name__ == '__main__':
    main()