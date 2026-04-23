#!/usr/bin/env python3
"""
Tree visualization tool for alifilter output
Generates an interactive HTML visualization of the Newick tree with:
- Branch-level coverage percentages displayed on branches
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


def analyze_branch_composition(node, coverage_data):
    """
    Analyze the composition of sequences under a node.
    Returns: (dominant_branch, num_sequences, coverage_percent, branch_breakdown)
    """
    leaves = collect_leaf_sequences(node)
    
    # Count branches
    branch_count = {}
    for leaf in leaves:
        branch = extract_branch_from_seq(leaf)
        branch_count[branch] = branch_count.get(branch, 0) + 1
    
    # Determine dominant branch
    dominant_branch = max(branch_count.items(), key=lambda x: x[1])[0] if branch_count else "Unknown"
    
    # Get coverage for dominant branch if available
    coverage_percent = 0.0
    if dominant_branch in coverage_data:
        coverage_percent = coverage_data[dominant_branch]['coverage_percent']
    
    return {
        'dominant_branch': dominant_branch,
        'sequence_count': len(leaves),
        'coverage_percent': coverage_percent,
        'branch_breakdown': branch_count,
        'leaves': leaves
    }


def generate_html_visualization(tree_path, coverage_json_path, output_html):
    """Generate interactive HTML visualization of the tree"""
    
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
    
    # Read tree
    with open(tree_path, 'r') as f:
        newick_str = f.read().strip()
    
    tree = parse_newick(newick_str)
    
    # Load coverage data (new format: {seq_name: {num, coverage_percent, removed_sequences}})
    with open(coverage_json_path, 'r') as f:
        sequence_coverage = json.load(f)
    
    # Build branch-level coverage statistics from sequence-level data
    branch_coverage = {}
    for seq_name, seq_data in sequence_coverage.items():
        branch = extract_branch_from_seq(seq_name)
        if branch not in branch_coverage:
            branch_coverage[branch] = {
                'total_num': 0,
                'total_coverage_percent': 0.0,
                'count': 0,
                'sequences': []
            }
        branch_coverage[branch]['total_num'] += seq_data['num']
        branch_coverage[branch]['total_coverage_percent'] += seq_data['coverage_percent']
        branch_coverage[branch]['count'] += 1
        branch_coverage[branch]['sequences'].append(seq_name)
    
    # Calculate average coverage per branch
    for branch in branch_coverage:
        branch_coverage[branch]['avg_coverage_percent'] = (
            branch_coverage[branch]['total_coverage_percent'] / 
            branch_coverage[branch]['count']
        ) if branch_coverage[branch]['count'] > 0 else 0.0
    
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
    
    # Generate SVG tree visualization with coverage on branches
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
                <circle cx="{x}" cy="{y}" r="5" fill="{color}" class="leaf-node {marker_class}" data-seq="{node.name}">
                    <title>{node.name}</title>
                </circle>
            ''')
            
            # Leaf label (shorter, cleaner)
            label = node.name.split('|')[-1] if '|' in node.name else node.name
            svg_elements.append(f'''
                <text x="{x + 10}" y="{y + 4}" font-size="11" class="leaf-label" data-seq="{node.name}">
                    {label}
                </text>
            ''')
            
            return svg_elements, x, y
        
        else:
            # Internal node - analyze composition
            analysis = analyze_branch_composition(node, branch_coverage)
            branch = analysis['dominant_branch']
            color = branch_colors.get(branch, '#CCCCCC')
            coverage_pct = analysis['coverage_percent']
            
            # Draw internal node
            svg_elements.append(f'''
                <circle cx="{x}" cy="{y}" r="4" fill="{color}" class="internal-node" opacity="0.7">
                    <title>{branch} ({analysis['sequence_count']} seqs, {coverage_pct:.1f}% coverage)</title>
                </circle>
            ''')
            
            # Process children
            num_children = len(node.children)
            child_y_spacing = dy / max(1, num_children)
            
            all_child_svg = []
            child_positions = []
            
            for i, child in enumerate(node.children):
                child_y = y + (i - (num_children - 1) / 2) * child_y_spacing
                child_positions.append((child_y, child))
            
            for child_y, child in child_positions:
                child_svg, child_x, child_y_ret = generate_tree_svg(
                    child, x + dx, child_y, dx, dy * 0.7, depth + 1
                )
                all_child_svg.extend(child_svg)
                
                # Draw branch line
                svg_elements.append(f'''
                    <line x1="{x}" y1="{y}" x2="{x + dx}" y2="{child_y}" 
                          stroke="{color}" stroke-width="2" opacity="0.8" class="branch-line"/>
                ''')
                
                # Add coverage percentage label on the branch (midpoint)
                mid_x = (x + x + dx) / 2
                mid_y = (y + child_y) / 2
                
                # Get coverage for the child's dominant branch
                child_analysis = analyze_branch_composition(child, branch_coverage)
                child_coverage = child_analysis['coverage_percent']
                
                if child_coverage > 0:
                    svg_elements.append(f'''
                        <text x="{mid_x}" y="{mid_y - 8}" 
                              font-size="12" font-weight="bold" text-anchor="middle" 
                              class="branch-coverage" fill="{color}">
                            {child_coverage:.1f}%
                        </text>
                    ''')
            
            svg_elements.extend(all_child_svg)
            
            return svg_elements, x, y
    
    svg_lines, _, _ = generate_tree_svg(tree, 80, 500, 180, 900)
    svg_content = '\n'.join(svg_lines)
    
    # Create legend
    legend_html = '''
        <div class="legend">
            <h3>Legend</h3>
            <div class="legend-section">
                <div class="legend-title">Branches</div>
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
                    <span class="legend-color" style="background-color: #96CEB4;"></span>
                    <span>SSL5</span>
                </div>
            </div>
            <div class="legend-section" style="margin-top: 16px; border-top: 1px solid #2a3a5e; padding-top: 12px;">
                <div class="legend-title">Sequence Types</div>
                <div class="legend-item">
                    <span class="legend-marker reference"></span>
                    <span>Always-Keep Reference</span>
                </div>
                <div class="legend-item">
                    <span class="legend-marker regular"></span>
                    <span>Regular Sequence</span>
                </div>
            </div>
        </div>
    '''
    
    # Create coverage table
    coverage_rows = ''
    for branch in sorted(sequence_coverage.keys()):
        stats = sequence_coverage[branch]
        num = stats['num']
        cov_pct = stats['coverage_percent']
        
        # Highlight reference sequences
        is_ref = branch in always_keep_set
        ref_marker = '⭐ ' if is_ref else ''
        
        coverage_rows += f'''
        <tr {'class="reference-row"' if is_ref else ''}>
            <td><strong>{ref_marker}{branch}</strong></td>
            <td class="numeric">{num}</td>
            <td class="numeric">{cov_pct:.2f}%</td>
        </tr>
        '''
    
    coverage_table = f'''
        <div class="coverage-section">
            <h3>Sequence Coverage Report</h3>
            <p style="font-size: 0.85em; color: #bdc3c7; margin-bottom: 12px;">
                Coverage metric: Number of sequences this representative accounts for, 
                divided by total sequences in the original alignment.
            </p>
            <table class="coverage-table">
                <thead>
                    <tr>
                        <th>Sequence</th>
                        <th>Represents</th>
                        <th>Coverage %</th>
                    </tr>
                </thead>
                <tbody>
                    {coverage_rows}
                </tbody>
            </table>
            
            <div class="coverage-explanation">
                <h4>📊 Understanding Coverage</h4>
                <div class="explanation-item">
                    <strong>Represents:</strong> The count of sequences this representative accounts for 
                    (1 = itself, plus any sequences merged with it due to &gt;90% identity)
                </div>
                <div class="explanation-item">
                    <strong>Coverage %:</strong> (Sequences Represented / Total Original Sequences) × 100
                </div>
                <div class="explanation-item" style="margin-top: 8px; padding-top: 8px; border-top: 1px solid #2a3a5e;">
                    <strong>Example:</strong> If a representative accounts for 25 sequences out of 1000 total, 
                    its coverage is 2.5%
                </div>
                <div class="explanation-item" style="color: #FFD700; margin-top: 8px;">
                    <strong>⭐ = Always-Keep Reference Sequence</strong> — Retained regardless of redundancy
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
            
            @import url('https://fonts.googleapis.com/css2?family=Playfair+Display:wght@700&family=Inter:wght@400;500;600&display=swap');
            
            body {{
                font-family: 'Inter', sans-serif;
                background: linear-gradient(135deg, #0a0e27 0%, #1a1f3a 100%);
                color: #ecf0f1;
                padding: 20px;
                min-height: 100vh;
            }}
            
            .container {{
                max-width: 1800px;
                margin: 0 auto;
                background: linear-gradient(135deg, #0f1629 0%, #1a2847 100%);
                border-radius: 16px;
                padding: 40px;
                box-shadow: 0 25px 80px rgba(0, 0, 0, 0.7);
                border: 1px solid rgba(0, 212, 255, 0.15);
            }}
            
            h1 {{
                text-align: center;
                margin-bottom: 8px;
                font-size: 2.8em;
                font-family: 'Playfair Display', serif;
                background: linear-gradient(135deg, #00d4ff 0%, #0099ff 100%);
                -webkit-background-clip: text;
                -webkit-text-fill-color: transparent;
                background-clip: text;
                letter-spacing: -1px;
            }}
            
            .subtitle {{
                text-align: center;
                color: #a0aec0;
                margin-bottom: 35px;
                font-size: 0.95em;
                font-weight: 500;
            }}
            
            .layout {{
                display: grid;
                grid-template-columns: 1fr 320px;
                gap: 30px;
                margin-bottom: 35px;
            }}
            
            .tree-container {{
                background: linear-gradient(135deg, #0d1b2a 0%, #1a2c42 100%);
                border-radius: 12px;
                overflow: auto;
                border: 2px solid rgba(0, 212, 255, 0.3);
                padding: 25px;
                max-height: 900px;
                box-shadow: inset 0 2px 8px rgba(0, 0, 0, 0.3);
            }}
            
            .tree-container svg {{
                display: block;
                margin: 0 auto;
                filter: drop-shadow(0 4px 12px rgba(0, 212, 255, 0.15));
            }}
            
            .leaf-node {{
                cursor: pointer;
                transition: all 0.3s ease;
                filter: drop-shadow(0 0 0px rgba(255, 255, 255, 0));
            }}
            
            .leaf-node:hover {{
                r: 7;
                filter: drop-shadow(0 0 10px currentColor);
            }}
            
            .reference-seq {{
                stroke: #FFD700;
                stroke-width: 2.5;
                filter: drop-shadow(0 0 8px rgba(255, 215, 0, 0.6));
            }}
            
            .reference-seq:hover {{
                stroke-width: 3;
                filter: drop-shadow(0 0 12px rgba(255, 215, 0, 0.8));
            }}
            
            .leaf-label {{
                fill: #d1d5db;
                cursor: pointer;
                transition: all 0.2s ease;
                font-size: 11px;
                font-weight: 500;
            }}
            
            .leaf-label:hover {{
                font-weight: 700;
                fill: #00d4ff;
                filter: drop-shadow(0 0 6px rgba(0, 212, 255, 0.5));
            }}
            
            .branch-line {{
                transition: all 0.3s ease;
            }}
            
            .branch-line:hover {{
                opacity: 1 !important;
                stroke-width: 2.5 !important;
                filter: drop-shadow(0 0 8px currentColor);
            }}
            
            .branch-coverage {{
                transition: all 0.3s ease;
                background: rgba(15, 52, 96, 0.9);
                padding: 2px 6px;
                border-radius: 3px;
                font-family: 'Inter', monospace;
                letter-spacing: 0.5px;
            }}
            
            .internal-node {{
                cursor: pointer;
                transition: all 0.2s ease;
                filter: drop-shadow(0 0 4px rgba(0, 0, 0, 0.3));
            }}
            
            .internal-node:hover {{
                r: 5;
                filter: drop-shadow(0 0 10px currentColor);
            }}
            
            .sidebar {{
                display: flex;
                flex-direction: column;
                gap: 24px;
            }}
            
            .legend {{
                background: linear-gradient(135deg, #0d1b2a 0%, #1a2c42 100%);
                border-radius: 12px;
                padding: 20px;
                border: 2px solid rgba(0, 212, 255, 0.3);
                box-shadow: inset 0 2px 8px rgba(0, 0, 0, 0.2);
            }}
            
            .legend h3 {{
                font-size: 1.25em;
                margin-bottom: 16px;
                color: #00d4ff;
                font-family: 'Playfair Display', serif;
                letter-spacing: -0.5px;
            }}
            
            .legend-section {{
                display: flex;
                flex-direction: column;
                gap: 10px;
            }}
            
            .legend-title {{
                font-size: 0.8em;
                text-transform: uppercase;
                letter-spacing: 1px;
                color: #7dd3fc;
                font-weight: 600;
                margin-bottom: 4px;
            }}
            
            .legend-item {{
                display: flex;
                align-items: center;
                gap: 10px;
                font-size: 0.9em;
                color: #cbd5e1;
                transition: all 0.2s ease;
                padding: 4px;
                border-radius: 4px;
            }}
            
            .legend-item:hover {{
                color: #00d4ff;
                background: rgba(0, 212, 255, 0.05);
            }}
            
            .legend-color {{
                width: 14px;
                height: 14px;
                border-radius: 3px;
                box-shadow: 0 2px 4px rgba(0, 0, 0, 0.2);
            }}
            
            .legend-marker {{
                width: 12px;
                height: 12px;
                border: 2px solid currentColor;
                border-radius: 50%;
                background: transparent;
            }}
            
            .legend-marker.reference {{
                border-color: #FFD700;
                box-shadow: 0 0 6px rgba(255, 215, 0, 0.4);
            }}
            
            .legend-marker.regular {{
                border-color: #94a3b8;
            }}
            
            .coverage-section {{
                background: linear-gradient(135deg, #0d1b2a 0%, #1a2c42 100%);
                border-radius: 12px;
                padding: 28px;
                border: 2px solid rgba(0, 212, 255, 0.3);
                grid-column: 1 / -1;
                box-shadow: inset 0 2px 8px rgba(0, 0, 0, 0.2);
            }}
            
            .coverage-section h3 {{
                font-size: 1.4em;
                margin-bottom: 12px;
                color: #00d4ff;
                font-family: 'Playfair Display', serif;
                letter-spacing: -0.5px;
            }}
            
            .coverage-section > p {{
                margin-bottom: 18px;
                line-height: 1.6;
            }}
            
            .coverage-table {{
                width: 100%;
                border-collapse: collapse;
                margin-bottom: 20px;
            }}
            
            .coverage-table thead {{
                background: rgba(0, 212, 255, 0.05);
                border-bottom: 2px solid rgba(0, 212, 255, 0.3);
            }}
            
            .coverage-table th {{
                padding: 12px 14px;
                text-align: left;
                font-weight: 600;
                color: #7dd3fc;
                font-size: 0.9em;
                letter-spacing: 0.5px;
                text-transform: uppercase;
            }}
            
            .coverage-table td {{
                padding: 11px 14px;
                border-bottom: 1px solid rgba(255, 255, 255, 0.05);
                color: #cbd5e1;
            }}
            
            .coverage-table td.numeric {{
                font-family: 'Inter', monospace;
                text-align: right;
            }}
            
            .coverage-table tr:hover {{
                background: rgba(0, 212, 255, 0.08);
            }}
            
            .coverage-table tbody tr.reference-row {{
                background: rgba(255, 215, 0, 0.08);
            }}
            
            .coverage-table tbody tr.reference-row:hover {{
                background: rgba(255, 215, 0, 0.15);
            }}
            
            .coverage-table tbody tr:last-child td {{
                border-bottom: none;
            }}
            
            .coverage-explanation {{
                font-size: 0.85em;
                color: #a0aec0;
                margin-top: 16px;
                padding: 16px;
                background: rgba(0, 212, 255, 0.05);
                border-left: 3px solid #00d4ff;
                border-radius: 6px;
            }}
            
            .coverage-explanation h4 {{
                color: #00d4ff;
                margin-bottom: 10px;
                font-size: 0.95em;
                font-weight: 600;
            }}
            
            .explanation-item {{
                padding: 6px 0;
                line-height: 1.5;
            }}
            
            @media (max-width: 1024px) {{
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
                
                .tree-container {{
                    max-height: 600px;
                }}
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>🧬 Phylogenetic Tree Analysis</h1>
            <div class="subtitle">Interactive visualization of sequence filtering and coverage metrics</div>
            
            <div class="layout">
                <div class="tree-container">
                    <svg viewBox="0 0 2000 1000" xmlns="http://www.w3.org/2000/svg">
                        <defs>
                            <style>
                                .leaf-node, .internal-node, .leaf-label, .branch-coverage {{ 
                                    transition: all 0.3s cubic-bezier(0.34, 1.56, 0.64, 1); 
                                }}
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
        
        <script>
            // Interactivity: highlight related sequences
            document.querySelectorAll('.leaf-node, .leaf-label').forEach(el => {{
                el.addEventListener('mouseenter', function() {{
                    const seq = this.getAttribute('data-seq');
                    if (seq) {{
                        document.querySelectorAll(`[data-seq="${{seq}}"]`).forEach(related => {{
                            related.style.opacity = '1';
                            related.style.filter = 'drop-shadow(0 0 12px rgba(0, 212, 255, 0.8))';
                        }});
                        document.querySelectorAll('[data-seq]').forEach(other => {{
                            if (other.getAttribute('data-seq') !== seq) {{
                                other.style.opacity = '0.3';
                            }}
                        }});
                    }}
                }});
                
                el.addEventListener('mouseleave', function() {{
                    document.querySelectorAll('[data-seq]').forEach(el => {{
                        el.style.opacity = '1';
                        el.style.filter = 'none';
                    }});
                }});
            }});
        </script>
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
        help='Path to sequence coverage JSON file (sequence_coverage.json)'
    )
    parser.add_argument(
        '-o', '--output', required=True,
        help='Output HTML file path'
    )
    
    args = parser.parse_args()
    
    # Generate visualization
    generate_html_visualization(args.tree, args.coverage, args.output)


if __name__ == '__main__':
    main()