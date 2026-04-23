#!/usr/bin/env python3
"""
Tree visualization tool for alifilter output
Generates an interactive HTML visualization of the Newick tree with:
- Branch-level coverage percentages
- Always-keep reference sequences highlighted
- Taxonomic coloring by branch
"""

import json
import os
import sys
import argparse
from pathlib import Path


def parse_newick(newick_str):
    """
    Parse a Newick format tree string.
    Returns a tree structure suitable for visualization.
    """
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
            # Internal node
            pos += 1  # skip '('
            while pos < len(newick_str) and newick_str[pos] != ')':
                if newick_str[pos] == ',':
                    pos += 1
                else:
                    child = parse_node()
                    node.children.append(child)
            
            pos += 1  # skip ')'
            
            # Parse node label (optional)
            label_start = pos
            while pos < len(newick_str) and newick_str[pos] not in '(),;:':
                pos += 1
            if pos > label_start:
                node.name = newick_str[label_start:pos]
        else:
            # Leaf node
            label_start = pos
            while pos < len(newick_str) and newick_str[pos] not in '(),;:':
                pos += 1
            node.name = newick_str[label_start:pos]
        
        # Parse distance (optional)
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


def extract_branch_from_seq(seq_name):
    """Extract branch (SSL3, SSL7, etc.) from sequence name"""
    if '|' in seq_name:
        return seq_name.split('|')[0]
    for prefix in ["SSL3", "SSL7", "SSL11", "SSL1", "SSL2"]:
        if seq_name.startswith(prefix):
            return prefix
    return "Unknown"


def collect_leaf_sequences(node):
    """Recursively collect all leaf sequences from a tree node"""
    if not node.children:
        return [node.name] if node.name else []
    
    leaves = []
    for child in node.children:
        leaves.extend(collect_leaf_sequences(child))
    return leaves


def analyze_branch_composition(node, branch_coverage):
    """
    Analyze the composition of sequences under a node.
    Returns: (branch_name, num_sequences, coverage_percent)
    """
    leaves = collect_leaf_sequences(node)
    
    # Count branches
    branch_count = {}
    for leaf in leaves:
        branch = extract_branch_from_seq(leaf)
        branch_count[branch] = branch_count.get(branch, 0) + 1
    
    # Determine dominant branch
    dominant_branch = max(branch_count.items(), key=lambda x: x[1])[0] if branch_count else "Unknown"
    
    return {
        'dominant_branch': dominant_branch,
        'sequence_count': len(leaves),
        'branch_breakdown': branch_count,
        'leaves': leaves
    }


def generate_html_visualization(tree_path, coverage_data, always_keep_data, output_html):
    """Generate interactive HTML visualization of the tree"""
    
    # Read tree
    with open(tree_path, 'r') as f:
        newick_str = f.read().strip()
    
    tree = parse_newick(newick_str)
    
    # Color mapping for branches
    branch_colors = {
        'SSL1': '#FF6B6B',
        'SSL2': '#4ECDC4',
        'SSL3': '#45B7D1',
        'SSL5': '#96CEB4',
        'SSL7': '#FFEAA7',
        'SSL11': '#DDA15E',
        'Unknown': '#CCCCCC'
    }
    
    # Prepare always-keep set
    always_keep_set = set(always_keep_data.get('found', []))
    
    # Generate SVG tree visualization
    def generate_tree_svg(node, x, y, dx, dy, depth=0):
        """Recursively generate SVG for tree nodes"""
        svg_elements = []
        
        if not node.children:
            # Leaf node
            branch = extract_branch_from_seq(node.name)
            color = branch_colors.get(branch, '#CCCCCC')
            
            is_reference = node.name in always_keep_set
            marker_class = 'reference-seq' if is_reference else 'regular-seq'
            
            # Leaf circle
            svg_elements.append(f'''
                <circle cx="{x}" cy="{y}" r="4" fill="{color}" class="leaf-node {marker_class}">
                    <title>{node.name}</title>
                </circle>
            ''')
            
            # Leaf label
            svg_elements.append(f'''
                <text x="{x + 8}" y="{y + 4}" font-size="10" class="leaf-label">
                    {node.name}
                </text>
            ''')
            
            return svg_elements, x, y
        
        else:
            # Internal node - analyze composition
            analysis = analyze_branch_composition(node, coverage_data)
            branch = analysis['dominant_branch']
            color = branch_colors.get(branch, '#CCCCCC')
            
            # Draw internal node
            svg_elements.append(f'''
                <circle cx="{x}" cy="{y}" r="3" fill="{color}" class="internal-node">
                    <title>{branch} ({analysis['sequence_count']} seqs)</title>
                </circle>
            ''')
            
            # Add coverage label if available
            if branch in coverage_data:
                coverage = coverage_data[branch]['coverage_percent']
                svg_elements.append(f'''
                    <text x="{x}" y="{y - 10}" font-size="9" font-weight="bold" 
                          text-anchor="middle" class="coverage-label">
                        {coverage:.1f}%
                    </text>
                ''')
            
            # Process children
            num_children = len(node.children)
            child_y_spacing = dy / max(1, num_children)
            
            child_ys = []
            all_child_svg = []
            
            for i, child in enumerate(node.children):
                child_y = y + (i - (num_children - 1) / 2) * child_y_spacing
                child_ys.append(child_y)
                
                child_svg, child_x, child_y = generate_tree_svg(
                    child, x + dx, child_y, dx, dy * 0.7, depth + 1
                )
                all_child_svg.extend(child_svg)
                
                # Draw branch line
                svg_elements.append(f'''
                    <line x1="{x}" y1="{y}" x2="{x + dx}" y2="{child_y}" 
                          stroke="{color}" stroke-width="1.5" opacity="0.7"/>
                ''')
            
            svg_elements.extend(all_child_svg)
            
            return svg_elements, x, y
    
    svg_lines, _, _ = generate_tree_svg(tree, 100, 500, 150, 800)
    svg_content = '\n'.join(svg_lines)
    
    # Create legend
    legend_html = '''
        <div class="legend">
            <h3>Legend</h3>
            <div class="legend-item">
                <span class="legend-color" style="background-color: #45B7D1;"></span>
                <span>SSL3</span>
            </div>
            <div class="legend-item">
                <span class="legend-color" style="background-color: #FFEAA7;"></span>
                <span>SSL7</span>
            </div>
            <div class="legend-item">
                <span class="legend-color" style="background-color: #DDA15E;"></span>
                <span>SSL11</span>
            </div>
            <div class="legend-item">
                <span class="legend-marker reference"></span>
                <span>Always-Keep Reference</span>
            </div>
        </div>
    '''
    
    # Create coverage table
    coverage_rows = ''
    for branch in sorted(coverage_data.keys()):
        stats = coverage_data[branch]
        num_reps = stats['num_representatives']
        total_removed = stats['total_removed_by_pid_to_reps']
        total_cluster = stats['total_cluster_size']
        coverage = stats['coverage_percent']
        
        coverage_rows += f'''
        <tr>
            <td><strong>{branch}</strong></td>
            <td>{num_reps}</td>
            <td>{total_removed}</td>
            <td>{total_cluster}</td>
            <td>{coverage:.2f}%</td>
        </tr>
        '''
    
    coverage_table = f'''
        <div class="coverage-section">
            <h3>Branch Coverage Statistics</h3>
            <p style="font-size: 0.9em; color: #bdc3c7; margin-bottom: 12px;">
                <strong>Coverage</strong> = Representatives / (Representatives + Sequences Removed by PID)
            </p>
            <table class="coverage-table">
                <thead>
                    <tr>
                        <th>Branch</th>
                        <th>Reps Kept</th>
                        <th>Removed by PID</th>
                        <th>Cluster Size</th>
                        <th>Coverage %</th>
                    </tr>
                </thead>
                <tbody>
                    {coverage_rows}
                </tbody>
            </table>
            
            <div class="rep-details">
                <h4>📊 How to Read This Table</h4>
                <div class="rep-item">
                    <strong>Reps Kept:</strong> Number of representative sequences retained for this branch
                </div>
                <div class="rep-item">
                    <strong>Removed by PID:</strong> Total sequences removed because they were &gt;90% identical to a kept representative
                </div>
                <div class="rep-item">
                    <strong>Cluster Size:</strong> Total sequences in clusters (reps + those removed by PID)
                </div>
                <div class="rep-item">
                    <strong>Coverage %:</strong> Reps / Cluster Size × 100 — represents the proportion of unique representatives
                </div>
                <div class="rep-item" style="margin-top: 8px; padding-top: 8px; border-top: 1px solid #2a3a5e;">
                    <strong>Example:</strong> If a branch has 100 reps and 200 removed by PID → 100 reps / (100+200) = 33.3% coverage
                </div>
            </div>
        </div>
    '''
    
    # Full HTML
    html_content = f'''
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Alifilter Tree Visualization</title>
        <style>
            * {{
                margin: 0;
                padding: 0;
                box-sizing: border-box;
            }}
            
            body {{
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                background: linear-gradient(135deg, #1a1a2e 0%, #16213e 100%);
                color: #ecf0f1;
                padding: 20px;
                min-height: 100vh;
            }}
            
            .container {{
                max-width: 1600px;
                margin: 0 auto;
                background: #0f3460;
                border-radius: 12px;
                padding: 30px;
                box-shadow: 0 20px 60px rgba(0, 0, 0, 0.5);
            }}
            
            h1 {{
                text-align: center;
                margin-bottom: 10px;
                font-size: 2.5em;
                background: linear-gradient(135deg, #00d4ff, #0099ff);
                -webkit-background-clip: text;
                -webkit-text-fill-color: transparent;
                background-clip: text;
            }}
            
            .subtitle {{
                text-align: center;
                color: #bdc3c7;
                margin-bottom: 30px;
                font-size: 0.95em;
            }}
            
            .layout {{
                display: grid;
                grid-template-columns: 1fr 300px;
                gap: 30px;
                margin-bottom: 30px;
            }}
            
            .tree-container {{
                background: #1a2942;
                border-radius: 8px;
                overflow: auto;
                border: 2px solid #00d4ff;
                padding: 20px;
                max-height: 800px;
            }}
            
            .tree-container svg {{
                display: block;
                margin: 0 auto;
                filter: drop-shadow(0 4px 6px rgba(0, 212, 255, 0.2));
            }}
            
            .leaf-node {{
                cursor: pointer;
                transition: r 0.2s ease;
            }}
            
            .leaf-node:hover {{
                r: 6;
                filter: drop-shadow(0 0 8px currentColor);
            }}
            
            .reference-seq {{
                stroke: #FFD700;
                stroke-width: 2;
                filter: drop-shadow(0 0 6px #FFD700);
            }}
            
            .leaf-label {{
                fill: #ecf0f1;
                cursor: pointer;
                transition: font-weight 0.2s ease;
            }}
            
            .leaf-label:hover {{
                font-weight: bold;
                fill: #00d4ff;
            }}
            
            .internal-node {{
                cursor: pointer;
                transition: r 0.2s ease;
                filter: drop-shadow(0 0 3px currentColor);
            }}
            
            .internal-node:hover {{
                r: 5;
            }}
            
            .coverage-label {{
                background: #0f3460;
                padding: 2px 4px;
                border-radius: 3px;
                fill: #FFD700;
                font-weight: bold;
                text-shadow: 0 0 4px rgba(0, 0, 0, 0.8);
            }}
            
            .sidebar {{
                display: flex;
                flex-direction: column;
                gap: 20px;
            }}
            
            .legend {{
                background: #1a2942;
                border-radius: 8px;
                padding: 15px;
                border: 2px solid #00d4ff;
            }}
            
            .legend h3 {{
                font-size: 1.2em;
                margin-bottom: 12px;
                color: #00d4ff;
            }}
            
            .legend-item {{
                display: flex;
                align-items: center;
                gap: 10px;
                margin-bottom: 8px;
                font-size: 0.9em;
            }}
            
            .legend-color {{
                width: 16px;
                height: 16px;
                border-radius: 3px;
            }}
            
            .legend-marker {{
                width: 12px;
                height: 12px;
                border: 2px solid #FFD700;
                border-radius: 50%;
            }}
            
            .coverage-section {{
                background: #1a2942;
                border-radius: 8px;
                padding: 20px;
                border: 2px solid #00d4ff;
                grid-column: 1 / -1;
            }}
            
            .coverage-section h3 {{
                font-size: 1.2em;
                margin-bottom: 15px;
                color: #00d4ff;
            }}
            
            .coverage-table {{
                width: 100%;
                border-collapse: collapse;
            }}
            
            .coverage-table thead {{
                background: #0f3460;
                border-bottom: 2px solid #00d4ff;
            }}
            
            .coverage-table th {{
                padding: 10px;
                text-align: left;
                font-weight: 600;
                color: #00d4ff;
            }}
            
            .coverage-table td {{
                padding: 8px 10px;
                border-bottom: 1px solid #2a3a5e;
            }}
            
            .coverage-table td {{
                padding: 8px 10px;
                border-bottom: 1px solid #2a3a5e;
            }}
            
            .coverage-table tr:hover {{
                background: #172445;
            }}
            
            .coverage-table tbody tr:last-child td {{
                border-bottom: none;
            }}
            
            .rep-details {{
                font-size: 0.85em;
                color: #bdc3c7;
                margin-top: 15px;
                padding: 15px;
                background: #0f3460;
                border-left: 3px solid #00d4ff;
                border-radius: 4px;
            }}
            
            .rep-details h4 {{
                color: #00d4ff;
                margin-bottom: 8px;
                font-size: 0.95em;
            }}
            
            .rep-item {{
                padding: 4px 0;
                line-height: 1.4;
            }}
                .layout {{
                    grid-template-columns: 1fr;
                }}
                
                .sidebar {{
                    grid-column: 1 / -1;
                    flex-direction: row;
                    gap: 20px;
                }}
                
                .legend {{
                    flex: 1;
                }}
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>🧬 Phylogenetic Tree Analysis</h1>
            <div class="subtitle">Interactive visualization of sequence filtering and coverage</div>
            
            <div class="layout">
                <div class="tree-container">
                    <svg viewBox="0 0 2000 1000" xmlns="http://www.w3.org/2000/svg">
                        <defs>
                            <style>
                                .leaf-node, .internal-node {{ transition: all 0.2s ease; }}
                            </style>
                        </defs>
                        {svg_content}
                    </svg>
                </div>
                
                <div class="sidebar">
                    {legend_html}
                </div>
            </div>
            
            {coverage_table}
        </div>
    </body>
    </html>
    '''
    
    with open(output_html, 'w') as f:
        f.write(html_content)
    
    print(f"✅ Visualization saved to: {output_html}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate interactive tree visualization from alifilter output'
    )
    parser.add_argument(
        '--tree', required=True,
        help='Path to Newick tree file (.tree)'
    )
    parser.add_argument(
        '--coverage', required=True,
        help='Path to branch coverage JSON file'
    )
    parser.add_argument(
        '--always-keep', required=True,
        help='Path to always-keep references JSON file'
    )
    parser.add_argument(
        '-o', '--output', required=True,
        help='Output HTML file path'
    )
    
    args = parser.parse_args()
    
    # Load data
    with open(args.coverage, 'r') as f:
        coverage_data = json.load(f)
    
    with open(args.always_keep, 'r') as f:
        always_keep_data = json.load(f)
    
    # Generate visualization
    generate_html_visualization(args.tree, coverage_data, always_keep_data, args.output)


if __name__ == '__main__':
    main()