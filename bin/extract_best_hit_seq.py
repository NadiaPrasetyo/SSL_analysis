#!/usr/bin/env python3
"""
parse_bath.py - Parse BATH (frameshift-aware HMM search) output files.

Extracts the best hit alignment for each query, translates the target DNA
codons to amino acids, and writes a FASTA file with the SSL query name
appended to the header.

Usage:
    python parse_bath.py <bath_output.txt> <output.fasta>
    python parse_bath.py <bath_output.txt>          # prints to stdout
"""

import sys
import re
import argparse
import glob
import os

# Standard codon table (table 1)
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def translate_codon(codon):
    codon = codon.upper().replace('U', 'T')
    return CODON_TABLE.get(codon, 'X')


def parse_bath_output(filepath):
    """
    Parse a BATH output file and return a list of best-hit records.

    Each record is a dict with:
        query       - query HMM name (e.g. Q2G2Y0|SSL7)
        target      - target sequence name
        evalue      - best hit E-value (float)
        score       - best hit score (float)
        ali_from    - alignment start on target (int)
        ali_to      - alignment end on target (int)
        aa_seq      - translated AA sequence from target codons
    """
    with open(filepath) as fh:
        content = fh.read()

    # Split into per-query blocks on "Query:" lines
    # Each block starts at a Query: line
    query_blocks = re.split(r'(?=^Query:)', content, flags=re.MULTILINE)

    records = []

    for block in query_blocks:
        if not block.strip():
            continue

        # ── Query name ────────────────────────────────────────────────────────
        query_match = re.search(r'^Query:\s+(\S+)', block, re.MULTILINE)
        if not query_match:
            continue
        query_name = query_match.group(1)

        # ── Scores table: grab the first (best) hit above inclusion threshold ─
        # Lines look like:
        #   3.3e-132  439.0  13.5  10000|654017|...  61022  60330  0  0
        score_pattern = re.compile(
            r'^\s+([\d.e+\-]+)\s+([\d.]+)\s+([\d.]+)\s+(\S+)\s+(\d+)\s+(\d+)',
            re.MULTILINE
        )
        score_matches = score_pattern.findall(block)

        # Filter out lines after "inclusion threshold" marker
        inclusion_idx = block.find('inclusion threshold')
        best_hit = None
        for m in score_matches:
            evalue_str, score_str, bias_str, target, start, end = m
            # check this match is before the inclusion threshold line
            match_pos = block.find(evalue_str + '  ' + score_str)
            if inclusion_idx != -1 and match_pos > inclusion_idx:
                continue
            best_hit = {
                'query': query_name,
                'target': target,
                'evalue': float(evalue_str),
                'score': float(score_str),
                'ali_from': int(start),
                'ali_to': int(end),
            }
            break  # first = best E-value

        if best_hit is None:
            continue

        # ── Alignment block for the best-scoring hit ─────────────────────────
        # Find the annotation section
        annot_section = re.search(
            r'Annotation for each hit.*', block, re.DOTALL
        )
        if not annot_section:
            continue
        annot_text = annot_section.group(0)

        # Within the annotation, find the hit block for the best target
        # Hit blocks start with ">> <target_name>"
        escaped_target = re.escape(best_hit['target'])
        hit_block_match = re.search(
            rf'^>>\s+{escaped_target}(.*?)(?=^>>|\Z)',
            annot_text,
            re.DOTALL | re.MULTILINE
        )
        if not hit_block_match:
            continue
        hit_block = hit_block_match.group(0)

        # Within one hit block there may be multiple sub-alignments (domains).
        # Grab only the first one (best score / lowest E-value).
        # Each sub-alignment has a score table row then an Alignment: section.
        # We split on "Alignment:" and take the first.
        alignment_parts = re.split(r'Alignment:', hit_block)
        if len(alignment_parts) < 2:
            continue
        first_alignment = alignment_parts[1]

        # ── Extract codons from the target line ───────────────────────────────
        # Target lines look like:
        #   10000|654017|...  61022  ATG AAA TTA ...  60969
        # They always contain the target name followed by a position number,
        # then space-separated uppercase 3-letter codons, then another number.
        codon_line_pattern = re.compile(
            rf'^\s+{escaped_target}\s+\d+\s+((?:[ACGT]{{3}}\s+)+)',
            re.MULTILINE
        )
        all_codons = []
        for line_match in codon_line_pattern.finditer(first_alignment):
            codons = line_match.group(1).split()
            all_codons.extend(codons)

        if not all_codons:
            print(f"[WARNING] No codons found for query {query_name}, "
                  f"target {best_hit['target']}", file=sys.stderr)
            continue

        aa_seq = ''.join(translate_codon(c) for c in all_codons)
        # Remove stop codons at the end (keep internal ones as X for visibility)
        aa_seq = aa_seq.rstrip('*')

        best_hit['aa_seq'] = aa_seq
        records.append(best_hit)

    return records


def write_fasta(records, out_fh):
    """Write records as FASTA to out_fh."""
    for rec in records:
        # Header: target sequence | query SSL name | coordinates | evalue
        header = (
            f">{rec['target']}"
            f"|{rec['query']} "
            f"ali_from={rec['ali_from']} ali_to={rec['ali_to']} "
            f"evalue={rec['evalue']}"
        )
        out_fh.write(header + '\n')
        # Wrap sequence at 60 chars
        seq = rec['aa_seq']
        for i in range(0, len(seq), 60):
            out_fh.write(seq[i:i+60] + '\n')


def main():
    parser = argparse.ArgumentParser(
        description="Extract best hit sequences from Bath output",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("-i", "--input", required=False, help="Bath output file")
    parser.add_argument("--input-dir", required=False, help="Directory containing Bath output files")
    parser.add_argument("-o", "--output", default="bath_best_hits.fasta", nargs='?', help="Output FASTA file (default: bath_best_hits.fasta)")
    parser.add_argument("--output-dir", default="./data", nargs='?', help="Output directory (default: data/)")

    args = parser.parse_args()

    if not args.input and not args.input_dir:
        parser.print_help()
        sys.exit(1)

    if args.input:
        records = parse_bath_output(args.input)
        
        if not records:
            print(f"[WARNING] No best hits extracted for {f}", file=sys.stderr)
            sys.exit(1)

        print(f"[INFO] Extracted {len(records)} best hit(s).", file=sys.stderr)
        output_file = args.output
        if output_file:
            with open(output_file, 'w') as fh:
                write_fasta(records, fh)
            print(f"[INFO] FASTA written to: {output_file}", file=sys.stderr)
        else:
            write_fasta(records, sys.stdout)
        sys.exit(0)

    if args.input_dir:
        input_files = glob.glob(os.path.join(args.input_dir, "*.out"))
        if not input_files:
            print("[WARNING] No input files found.", file=sys.stderr)
            sys.exit(1)
        records = []
        for f in input_files:
            records = (parse_bath_output(f))
            
            print(f"[INFO] Extracted {len(records)} best hit(s).", file=sys.stderr)

            f_stem = f.split("/")[-1].split("_")[0]
            os.makedirs(args.output_dir, exist_ok=True)
            output_file = os.path.join(args.output_dir, f"{f_stem}_best_hits.fasta")
            with open(output_file, 'w') as fh:
                write_fasta(records, fh)
            print(f"[INFO] FASTA written to: {output_file}", file=sys.stderr)
        sys.exit(0)


if __name__ == '__main__':
    main()
