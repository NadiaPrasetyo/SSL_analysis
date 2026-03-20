#!/usr/bin/env python3
"""
annotate_alignment.py

Parses raw NetMHCpan/NetMHCIIpan output files, selects the top N peptides
by score (lowest %Rank), and annotates a ClustalW alignment with colour-coded
peptide boxes rendered as an HTML file.

Usage:
    python annotate_alignment.py \
        --alignment alignment.clustal \
        --mhc-dir  /path/to/mhc_files \
        --top-n    10 \
        --output   annotated_alignment.html \
        [--mhcii]                   # flag if files are MHC-II
"""

import argparse
import os
import re
import statistics
from collections import defaultdict
from pathlib import Path


# ─────────────────────────────────────────────────────────────────────────────
# 1.  MHC file parser  (adapted from parse_mhc_dir in the original codebase)
# ─────────────────────────────────────────────────────────────────────────────

def parse_mhc_files(directory: str, is_mhcii: bool = False):
    """
    Parse raw NetMHCpan / NetMHCIIpan prediction files.

    Returns a flat list of peptide records:
        {
          "accession": str,   # e.g. "SSL3|CC1"
          "peptide":   str,   # amino-acid sequence
          "pos":       int,   # 1-based position in the source sequence
          "score":     float, # EL-score
          "rank":      float, # %Rank
          "binding":   str,   # "SB", "WB", or ""
          "allele":    str,
        }
    """
    prefix = "mhcii" if is_mhcii else "mhci"
    records = []

    try:
        files = [f for f in os.listdir(directory) if not f.endswith(".xls")]
    except Exception as e:
        raise SystemExit(f"Cannot list directory '{directory}': {e}")

    for fname in files:
        path = os.path.join(directory, fname)
        try:
            with open(path) as fh:
                lines = fh.readlines()
        except Exception as e:
            print(f"  [warn] Cannot read {fname}: {e}")
            continue

        for i, line in enumerate(lines):
            # Skip header / separator / empty lines
            if (line.startswith("#") or line.startswith("-") or
                    not line.strip() or "Pos" in line):
                continue
            # Only data lines contain a numeric first column
            parts = line.split()
            if len(parts) < 10:
                continue
            try:
                pos = int(parts[0])
            except ValueError:
                continue

            try:
                if prefix == "mhci":
                    # Col indices for NetMHCpan tabular output
                    allele   = parts[1]
                    peptide  = parts[2]
                    identity = parts[10]          # e.g. SSL3_CC1
                    score    = float(parts[11])
                    rank     = float(parts[12])
                    binding  = parts[14] if len(parts) > 14 else ""
                else:
                    # NetMHCIIpan
                    allele   = parts[1]
                    peptide  = parts[2]
                    identity = parts[7]
                    score    = float(parts[8])
                    rank     = float(parts[9])
                    binding  = parts[12] if len(parts) > 12 else ""

                # Normalise identity  "SSL3_CC1" → "SSL3|CC1"
                accession = identity.replace("_", "|", 1)

                records.append({
                    "accession": accession,
                    "peptide":   peptide,
                    "pos":       pos,          # 1-based
                    "score":     score,
                    "rank":      rank,
                    "binding":   binding,
                    "allele":    allele,
                })
            except (ValueError, IndexError) as e:
                print(f"  [debug] Skipping line {i} in {fname}: {e}")

    return records


# ─────────────────────────────────────────────────────────────────────────────
# 2.  Select top-N peptides per accession  (lowest %Rank = best binders)
# ─────────────────────────────────────────────────────────────────────────────

def select_top_n(records, n: int):
    """
    For every accession, pick the N peptides with the lowest %Rank.
    Returns a dict:  accession → [record, ...]
    """
    by_accession = defaultdict(list)
    for r in records:
        by_accession[r["accession"]].append(r)

    top = {}
    for acc, recs in by_accession.items():
        sorted_recs = sorted(recs, key=lambda x: x["rank"])
        top[acc] = sorted_recs[:n]
    return top


# ─────────────────────────────────────────────────────────────────────────────
# 3.  Clustal alignment parser
# ─────────────────────────────────────────────────────────────────────────────

def parse_clustal(path: str):
    """
    Returns:
        sequences : dict  name → full ungapped sequence string (no spaces)
        blocks    : list of (name_col_width, [(name, seq_chunk), ...], conservation_line)
        order     : list of names in appearance order
    """
    with open(path) as fh:
        raw = fh.read()

    lines = raw.splitlines()
    # Skip the header line
    data_lines = [l for l in lines[1:] if l.strip()]

    # Group into blocks separated by the conservation line (starts with space or *)
    blocks = []          # list of lists-of-(name, chunk)
    current = []
    conservation_lines = []

    for line in data_lines:
        # Conservation / ruler line: starts with spaces but no sequence identifier
        if re.match(r'^\s+[\s\*\:\.\-]*$', line):
            if current:
                blocks.append(current)
                conservation_lines.append(line)
                current = []
            continue
        parts = line.split()
        if len(parts) >= 2:
            name  = parts[0]
            chunk = parts[1]
            current.append((name, chunk))

    if current:
        blocks.append(current)
        conservation_lines.append("")

    # Build full sequences (with gaps)
    order = []
    sequences = {}
    for name, chunk in blocks[0]:
        order.append(name)
        sequences[name] = ""

    for block in blocks:
        for name, chunk in block:
            sequences[name] += chunk

    return sequences, blocks, conservation_lines, order


# ─────────────────────────────────────────────────────────────────────────────
# 4.  Map peptide positions onto the gapped alignment
# ─────────────────────────────────────────────────────────────────────────────

def ungapped_to_gapped_index(gapped_seq: str, ungapped_pos: int) -> int:
    """
    Given a 1-based position in the ungapped sequence, return the 0-based
    index in the gapped (alignment) sequence.
    """
    count = 0
    for i, ch in enumerate(gapped_seq):
        if ch != '-':
            count += 1
            if count == ungapped_pos:
                return i
    return len(gapped_seq)  # past end


def build_annotation_map(top_peptides, sequences):
    """
    Returns a dict:
        (accession, gapped_start, gapped_end) → peptide_record

    gapped_start / gapped_end are 0-based indices into the full gapped sequence.
    """
    annotations = []
    for acc, recs in top_peptides.items():
        gapped = sequences.get(acc)
        if gapped is None:
            print(f"  [warn] Accession '{acc}' not found in alignment — skipping")
            continue
        for rec in recs:
            start_g = ungapped_to_gapped_index(gapped, rec["pos"])
            end_g   = ungapped_to_gapped_index(gapped, rec["pos"] + len(rec["peptide"]) - 1)
            annotations.append({
                "accession":    acc,
                "gapped_start": start_g,
                "gapped_end":   end_g,
                "record":       rec,
            })
    return annotations


# ─────────────────────────────────────────────────────────────────────────────
# 5.  HTML renderer
# ─────────────────────────────────────────────────────────────────────────────

def make_family_colour_map(accessions):
    """
    Groups accessions by their protein-family prefix (the part before '|', e.g. 'SSL3').
    Each family gets a distinct hue band; within a family, shades vary by lightness
    so every individual sequence has a clearly distinguishable but obviously related colour.

    Returns:
        colour_map : dict  accession -> CSS colour string
        family_hues: dict  family_prefix -> hue (int, 0-359)
    """
    from collections import defaultdict

    # Collect families while preserving appearance order
    family_of = {}
    families_seen = []
    for acc in sorted(accessions):
        prefix = acc.split("|")[0] if "|" in acc else acc.split("_")[0]
        family_of[acc] = prefix
        if prefix not in families_seen:
            families_seen.append(prefix)

    n_families = len(families_seen)

    # Spread family hues evenly, starting at 210 (blue) for a pleasant default
    family_hues = {}
    for i, fam in enumerate(families_seen):
        hue = int((210 + 360 * i / n_families) % 360)
        family_hues[fam] = hue

    # Within each family vary lightness: darker members are more saturated
    family_members = defaultdict(list)
    for acc in sorted(accessions):
        family_members[family_of[acc]].append(acc)

    colour_map = {}
    for fam, members in family_members.items():
        hue = family_hues[fam]
        n   = len(members)
        for j, acc in enumerate(members):
            # Lightness: spread from 35% (rich) to 60% (pastel) across members
            lightness  = int(35 + 25 * j / max(n - 1, 1))
            saturation = 80 if lightness < 50 else 65
            colour_map[acc] = f"hsl({hue},{saturation}%,{lightness}%)"

    return colour_map, family_hues


def render_html(sequences, blocks, conservation_lines, order,
                annotations, top_peptides, output_path):

    # Assign family-grouped colours to each accession that has annotations
    annotated_accs = sorted({a["accession"] for a in annotations})
    colour_map, family_hues = make_family_colour_map(annotated_accs)

    # Build a per-(accession, col) lookup for quick highlight checks
    # Key: (accession, col_index) → annotation record
    highlight = {}
    for ann in annotations:
        acc = ann["accession"]
        for col in range(ann["gapped_start"], ann["gapped_end"] + 1):
            key = (acc, col)
            # Keep the best-ranked peptide if overlapping
            if key not in highlight or ann["record"]["rank"] < highlight[key]["record"]["rank"]:
                highlight[key] = ann

    # ── Build the legend HTML (grouped by family) ────────────────────────
    from collections import defaultdict as _dd
    family_accs = _dd(list)
    for acc in annotated_accs:
        prefix = acc.split("|")[0] if "|" in acc else acc.split("_")[0]
        family_accs[prefix].append(acc)

    legend_rows = []
    for family, members in family_accs.items():
        fam_hue    = family_hues.get(family, 210)
        fam_colour = f"hsl({fam_hue},80%,42%)"
        member_blocks = []
        for acc in members:
            colour = colour_map[acc]
            recs   = top_peptides.get(acc, [])
            rows   = "".join(
                f"<tr><td>{r['peptide']}</td><td>{r['allele']}</td>"
                f"<td>{r['rank']:.3f}</td><td>{r['binding'] or '—'}</td></tr>"
                for r in recs
            )
            member_blocks.append(f"""
            <div class="legend-member">
              <div class="legend-swatch" style="background:{colour}"></div>
              <strong>{acc}</strong>
              <table>
                <thead><tr><th>Peptide</th><th>Allele</th><th>%Rank</th><th>Level</th></tr></thead>
                <tbody>{rows}</tbody>
              </table>
            </div>""")
        legend_rows.append(f"""
        <div class="legend-family" style="border-left:4px solid {fam_colour}">
          <div class="legend-family-header" style="color:{fam_colour}">{family}</div>
          {"".join(member_blocks)}
        </div>""")

    # ── Build alignment HTML ───────────────────────────────────────────────
    block_htmls = []
    col_offset  = 0   # running alignment column offset

    for b_idx, block in enumerate(blocks):
        cons_line = conservation_lines[b_idx] if b_idx < len(conservation_lines) else ""

        # Determine the width of the name column
        name_width = max(len(name) for name, _ in block) + 2

        seq_lines_html = []
        block_len = len(block[0][1]) if block else 0

        for name, chunk in block:
            # Build each character with possible highlight
            chars_html = []
            for local_i, ch in enumerate(chunk):
                col = col_offset + local_i
                key = (name, col)
                if key in highlight:
                    ann    = highlight[key]
                    colour = colour_map[ann["accession"]]
                    tip    = (f"Peptide: {ann['record']['peptide']} | "
                              f"Allele: {ann['record']['allele']} | "
                              f"%Rank: {ann['record']['rank']:.3f} | "
                              f"Level: {ann['record']['binding'] or 'NB'}")
                    chars_html.append(
                        f'<span class="pep" style="background:{colour}" title="{tip}">{ch}</span>'
                    )
                else:
                    chars_html.append(f'<span class="aa">{ch}</span>')

            name_pad = name.ljust(name_width)
            seq_html  = "".join(chars_html)
            seq_lines_html.append(
                f'<div class="seq-row"><span class="seq-name">{name_pad}</span>'
                f'<span class="seq-body">{seq_html}</span></div>'
            )

        # Conservation line
        cons_html = (f'<div class="cons-row">'
                     f'<span class="seq-name">{"".ljust(name_width)}</span>'
                     f'<span class="seq-body cons">{cons_line.strip()}</span>'
                     f'</div>')

        block_htmls.append(
            f'<div class="aln-block">{"".join(seq_lines_html)}{cons_html}</div>'
        )
        col_offset += block_len

    # ── Assemble full page ─────────────────────────────────────────────────
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>MHC Peptide Alignment Annotation</title>
<style>
  :root {{
    --bg:      #f5f7fa;
    --surface: #ffffff;
    --border:  #dde1ea;
    --text:    #1a202c;
    --muted:   #64748b;
    --mono:    'Courier New', monospace;
  }}
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{
    background: var(--bg);
    color: var(--text);
    font-family: 'Segoe UI', system-ui, sans-serif;
    padding: 2rem;
  }}
  h1 {{ font-size: 1.6rem; margin-bottom: 0.25rem; color: #1a202c; }}
  .subtitle {{ color: var(--muted); font-size: 0.9rem; margin-bottom: 2rem; }}

  /* ── Legend ── */
  h2 {{ font-size: 1rem; text-transform: uppercase;
        letter-spacing: 0.08em; color: var(--muted);
        margin-bottom: 1rem; }}
  .legend-grid {{
    display: flex; flex-wrap: wrap; gap: 1.4rem; margin-bottom: 2.5rem;
  }}
  /* Family container */
  .legend-family {{
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 0.9rem 1.1rem;
    min-width: 300px;
    flex: 1;
    display: flex;
    flex-direction: column;
    gap: 0.7rem;
  }}
  .legend-family-header {{
    font-size: 0.7rem;
    font-weight: 800;
    text-transform: uppercase;
    letter-spacing: 0.12em;
    margin-bottom: 0.1rem;
  }}
  /* Per-accession block inside a family */
  .legend-member {{
    background: var(--bg);
    border: 1px solid var(--border);
    border-radius: 6px;
    padding: 0.55rem 0.75rem;
  }}
  .legend-swatch {{
    display: inline-block; width: 11px; height: 11px;
    border-radius: 3px; margin-right: 6px; vertical-align: middle;
  }}
  .legend-member table {{
    width: 100%; border-collapse: collapse;
    font-size: 0.76rem; margin-top: 0.45rem;
  }}
  .legend-member th {{
    color: var(--muted); text-align: left;
    padding: 2px 6px; border-bottom: 1px solid var(--border);
  }}
  .legend-member td {{ padding: 2px 6px; font-family: var(--mono); }}
  .legend-member tr:hover td {{ background: rgba(0,0,0,0.04); }}

  /* ── Alignment ── */
  .alignment-wrap {{
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 1.5rem;
    overflow-x: auto;
  }}
  .aln-block {{ margin-bottom: 1.2rem; }}
  .seq-row, .cons-row {{
    display: flex; align-items: baseline;
    line-height: 1.55;
  }}
  .seq-name {{
    font-family: var(--mono);
    font-size: 0.82rem;
    color: var(--muted);
    white-space: pre;
    flex-shrink: 0;
    padding-right: 0.5rem;
  }}
  .seq-body {{
    font-family: var(--mono);
    font-size: 0.88rem;
    letter-spacing: 0.04em;
    white-space: pre;
  }}
  .aa {{ color: var(--text); }}
  .pep {{
    border-radius: 2px;
    color: #fff;
    font-weight: 700;
    cursor: default;
    transition: opacity 0.15s;
  }}
  .pep:hover {{ opacity: 0.75; }}
  .cons {{ color: #2563eb; font-weight: 500; }}
</style>
</head>
<body>

<h1>MHC Epitope Alignment Annotation</h1>
<p class="subtitle">Top-ranked predicted binders highlighted on the multiple sequence alignment.</p>

<h2>Legend — Top Peptides per Sequence</h2>
<div class="legend-grid">
{"".join(legend_rows)}
</div>

<h2>Annotated Alignment</h2>
<div class="alignment-wrap">
{"".join(block_htmls)}
</div>

</body>
</html>
"""

    with open(output_path, "w") as fh:
        fh.write(html)
    print(f"\n✓ Annotated alignment written to: {output_path}")


# ─────────────────────────────────────────────────────────────────────────────
# 6.  CLI entry-point
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Annotate a ClustalW alignment with top-N MHC peptide binders."
    )
    parser.add_argument("--alignment", required=True,
                        help="Path to the ClustalW alignment file.")
    parser.add_argument("--mhc-dir", required=True,
                        help="Directory containing raw NetMHCpan/NetMHCIIpan output files.")
    parser.add_argument("--top-n", type=int, default=10,
                        help="Number of top peptides (by %%Rank) to annotate per sequence. Default: 10")
    parser.add_argument("--output", default="annotated_alignment.html",
                        help="Output HTML file path. Default: annotated_alignment.html")
    parser.add_argument("--mhcii", action="store_true",
                        help="Set this flag if the MHC files are MHC-II predictions (NetMHCIIpan format).")
    args = parser.parse_args()

    # ── Parse MHC predictions ──────────────────────────────────────────────
    print(f"[1/4] Parsing MHC prediction files in '{args.mhc_dir}' …")
    records = parse_mhc_files(args.mhc_dir, is_mhcii=args.mhcii)
    print(f"      → {len(records)} peptide records loaded.")

    if not records:
        raise SystemExit("No valid peptide records found. Check your MHC directory and file format.")

    # ── Select top-N per accession ─────────────────────────────────────────
    print(f"[2/4] Selecting top {args.top_n} peptides per accession …")
    top_peptides = select_top_n(records, args.top_n)
    for acc, recs in top_peptides.items():
        print(f"      {acc}: {len(recs)} peptides selected "
              f"(best %Rank = {recs[0]['rank']:.3f})")

    # ── Parse alignment ────────────────────────────────────────────────────
    print(f"[3/4] Parsing alignment '{args.alignment}' …")
    sequences, blocks, conservation_lines, order = parse_clustal(args.alignment)
    print(f"      → {len(sequences)} sequences, {len(blocks)} alignment blocks.")

    # ── Map peptides → alignment columns ──────────────────────────────────
    print("[4/4] Mapping peptides onto alignment columns …")
    annotations = build_annotation_map(top_peptides, sequences)
    print(f"      → {len(annotations)} peptide–column mappings built.")

    # ── Render HTML ────────────────────────────────────────────────────────
    render_html(sequences, blocks, conservation_lines, order,
                annotations, top_peptides, args.output)


if __name__ == "__main__":
    main()