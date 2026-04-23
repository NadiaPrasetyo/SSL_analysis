#!/usr/bin/env python3
# Usage: esl_alifilter.py <maxid> <msafile>
# e.g.:  esl_alifilter.py 0.90 my_alignment.sto
import sys
import subprocess
import tempfile
import os
import argparse
import logging
import shutil
import json
from collections import defaultdict

def setup_logging(verbose, output_dir):
    log_level = logging.DEBUG if verbose else logging.INFO
    log_file = os.path.join(output_dir, "alifilter.log")
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file, mode='w') if verbose else logging.NullHandler(),
            logging.StreamHandler(sys.stdout)
        ]
    )
    if verbose:
        logging.info("Logging initialized. Log file: %s", log_file)

def get_year(seq_name):
    """Extract year from header format: accession_country_year_ST_protein"""
    parts = seq_name.split("_")
    try:
        return int(parts[2])
    except (IndexError, ValueError):
        return 0  # fallback if parsing fails

def get_branch_from_sequence(seq_name):
    """Extract branch identifier from sequence name (e.g., SSL3, SSL7, SSL11)"""
    if '|' in seq_name:
        branch = seq_name.split('|')[0]
        return branch
    # fallback: try to extract from beginning
    for prefix in ["SSL3", "SSL7", "SSL11", "SSL1", "SSL2"]:
        if seq_name.startswith(prefix):
            return prefix
    return "Unknown"

def get_seqs_in_stockholm(stockholm_file):
    """Extract all sequence names from a Stockholm format file.
    Returns both full names (with prefix) and base names (without prefix)."""
    full_names = set()
    base_names = {}  # maps base name -> full name (with prefix)
    
    with open(stockholm_file, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            # In Stockholm format, sequences are like: "NAME  SEQUENCE"
            # Skip comments and blank lines
            if not line or line.startswith('#') or line.startswith('//'):
                continue
            # Split on whitespace
            parts = line.split()
            if len(parts) >= 1 and not line.startswith(' '):
                # This is a sequence line
                full_name = parts[0]
                full_names.add(full_name)
                
                # Extract base name (remove numeric prefix like "000|")
                # Format: "000|original_name"
                if '|' in full_name:
                    base_name = '|'.join(full_name.split('|')[1:])
                else:
                    base_name = full_name
                    
                base_names[base_name] = full_name
    
    return full_names, base_names

def parse_newick_tree(tree_file):
    """Parse a Newick format tree file and extract branch structure"""
    with open(tree_file, 'r') as f:
        newick_str = f.read().strip()
    
    # Simple Newick parser to extract branch names and structure
    branches = {}
    
    # Extract leaf node labels (sequences)
    import re
    # Match patterns like "name:distance" or just "name"
    pattern = r'([A-Za-z0-9|_\-\.]+)(?::[\d\.e\-]+)?'
    matches = re.findall(pattern, newick_str)
    
    for match in matches:
        if match and not match.startswith('(') and not match.startswith(')'):
            branch = get_branch_from_sequence(match)
            if branch not in branches:
                branches[branch] = {'total': 0, 'kept': 0, 'removed': 0, 'sequences': []}
            branches[branch]['sequences'].append(match)
    
    return branches

def main(tool_root, maxid, msafile, output_file, verbose=False):
    alipid_path      = os.path.join(tool_root, "easel", "miniapps", "esl-alipid")
    alimanip_path    = os.path.join(tool_root, "easel", "miniapps", "esl-alimanip")
    reformat_path = os.path.join(tool_root, "easel", "miniapps", "esl-reformat")
    output_dir = os.path.dirname(output_file)
    if not output_dir:
        output_dir = "."

    # Step 1: convert input AFA -> Stockholm first (to get canonical sequence names)
    os.makedirs(os.path.join(output_dir, "temp"), exist_ok=True)
    tmp_sto = os.path.join(output_dir, "temp", "tmp.sto")
    if verbose:
        logging.info(f"Converting {msafile} to Stockholm format...")
    with open(tmp_sto, "w") as out:
        subprocess.run(
            [str(reformat_path), "--informat", "afa", "stockholm", msafile],
            stdout=out, check=True
        )

    # Get all sequences actually in the Stockholm file
    full_names, base_names_map = get_seqs_in_stockholm(tmp_sto)
    if verbose:
        logging.info(f"Found {len(full_names)} sequences in Stockholm alignment")
        logging.debug(f"Sample full names: {list(full_names)[:5]}")
        logging.debug(f"Sample base name mappings: {dict(list(base_names_map.items())[:5])}")

    # Track sequences by branch for coverage calculation
    branch_stats = defaultdict(lambda: {'total': 0, 'kept': [], 'removed': []})
    
    for full_name in full_names:
        branch = get_branch_from_sequence(full_name)
        branch_stats[branch]['total'] += 1

    # Run esl-alipid to get all pairwise IDs
    if verbose:
        logging.info("Running esl-alipid to get pairwise identities...")
    result = subprocess.run(
        [str(alipid_path), "--amino", "--informat", "afa", msafile],
        capture_output=True, text=True, check=True
    )
    
    # Parse pairwise IDs; greedily remove the "worse" of any pair above threshold
    # with preference given to newer sequences and against unknown sequences
    
    # Track which sequences had high PID pairs with which reps
    pid_pairs = defaultdict(list)  # Maps rep_seq -> [(high_pid_seq, pid), ...]
    
    removed = set()
    all_seqs_ordered = []
    all_seqs_seen = set()
    skipped_count = 0
    processed_count = 0
    for line in result.stdout.splitlines():
        if line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 3:
            continue
        seq1, seq2, pid = fields[0], fields[1], float(fields[2])
        
        # Convert base names to full names (with prefix)
        full_seq1 = base_names_map.get(seq1, seq1)
        full_seq2 = base_names_map.get(seq2, seq2)
        
        # Only process if both sequences exist in the alignment
        if full_seq1 not in full_names or full_seq2 not in full_names:
            if verbose:
                logging.debug(f"Skipping pair: {seq1} -> {full_seq1} vs {seq2} -> {full_seq2} (one or both not in alignment)")
            skipped_count += 1
            continue
        
        processed_count += 1
        
        for s in (full_seq1, full_seq2):
            if s not in all_seqs_seen:
                all_seqs_ordered.append(s)
                all_seqs_seen.add(s)
        
        if pid >= maxid * 100.0:
            year1 = get_year(seq1)
            year2 = get_year(seq2)
            has_unknown_1 = "unknown" in full_seq1.lower()
            has_unknown_2 = "unknown" in full_seq2.lower()
            
            # Decision logic: prefer sequences without "unknown"
            if has_unknown_1 and not has_unknown_2:
                to_remove = full_seq1  # prefer seq2 (no unknown)
                to_keep = full_seq2
            elif has_unknown_2 and not has_unknown_1:
                to_remove = full_seq2  # prefer seq1 (no unknown)
                to_keep = full_seq1
            elif year1 > 2016 and year2 <= 2016:
                to_remove = full_seq2  # prefer seq1 (newer)
                to_keep = full_seq1
            elif year2 > 2016 and year1 <= 2016:
                to_remove = full_seq1  # prefer seq2 (newer)
                to_keep = full_seq2
            else:
                to_remove = full_seq2  # both same era: fall back to original behaviour
                to_keep = full_seq1
            
            if to_remove not in removed:
                if to_keep not in removed:
                    removed.add(to_remove)
                    # Track which sequences are removed relative to this kept sequence
                    pid_pairs[to_keep].append((to_remove, pid))
                    # Track removal by branch
                    branch = get_branch_from_sequence(to_remove)
                    branch_stats[branch]['removed'].append(to_remove)
                    if verbose:
                        logging.debug(f"Marking for removal: {to_remove} (PID: {pid:.1f}% with {to_keep})")

    if verbose:
        logging.info(f"Processed {processed_count} pairs where sequences are in alignment")
        logging.info(f"Skipped {skipped_count} pairs where sequences weren't in alignment")
        logging.info(f"Removed {len(removed)} sequences based on pairwise identity")
    
    # Track kept sequences by branch
    keep = [s for s in all_seqs_ordered if s not in removed]
    for seq in keep:
        branch = get_branch_from_sequence(seq)
        branch_stats[branch]['kept'].append(seq)

    always_keep_SSL3 = [
        "SSL3|CC1", "SSL3|CC5", "SSL3|CC8", "SSL3|CC22", "SSL3|CC30", "SSL3|CC93"
    ]
    always_keep_SSL7 = [
        "SSL7|CC1", "SSL7|CC5", "SSL7|CC8", "SSL7|CC22", "SSL7|CC30", "SSL7|CC93"
    ]
    always_keep_SSL11 = [
        "SSL11|CC1", "SSL11|CC5", "SSL11|CC8", "SSL11|CC22", "SSL11|CC30", "SSL11|CC93"
    ]

    # Track which always-keep sequences are kept
    always_keep_found = set()

    # Add specific sequences to keep list based on output file name
    # ONLY if they actually exist in the alignment
    added_count = 0
    if "SSL3" in output_file:
        for s in always_keep_SSL3:
            if s in base_names_map:
                full_s = base_names_map[s]
                if full_s not in keep:
                    keep.append(full_s)
                    added_count += 1
                always_keep_found.add(s)
    elif "SSL7" in output_file:
        for s in always_keep_SSL7:
            if s in base_names_map:
                full_s = base_names_map[s]
                if full_s not in keep:
                    keep.append(full_s)
                    added_count += 1
                always_keep_found.add(s)
    elif "SSL11" in output_file:
        for s in always_keep_SSL11:
            if s in base_names_map:
                full_s = base_names_map[s]
                if full_s not in keep:
                    keep.append(full_s)
                    added_count += 1
                always_keep_found.add(s)
    
    if verbose:
        logging.info(f"Added {added_count} reference sequences to keep list")

    # Deduplicate and clean
    keep = list(set(keep))
    keep = sorted([s.strip() for s in keep if s and s.strip()])
    
    if verbose:
        logging.info(f"Final keep list: {len(keep)} sequences")
        if len(keep) == 0:
            logging.error("ERROR: Keep list is empty! No sequences will be retained.")
            logging.error(f"Alignment has {len(full_names)} sequences")
            return
        logging.info(f"Sample keep sequences: {keep[:5]}")

    # Write keep-list to a persistent intermediate file (not temp)
    to_keep_file = os.path.join(output_dir, "temp", "to_keep.txt")
    
    # Write carefully - one sequence per line, no extra whitespace
    with open(to_keep_file, "w") as f:
        for seq in keep:
            # Ensure no trailing whitespace
            f.write(seq.strip() + "\n")
    
    if verbose:
        logging.info(f"Wrote keep list to {to_keep_file}")
        # Verify the file was written correctly
        with open(to_keep_file, "r") as f:
            file_contents = f.read()
            written_lines = [line for line in file_contents.split('\n') if line.strip()]
        logging.info(f"Verified {len(written_lines)} non-empty sequences in keep list file")
        
        # Show file size and content sample
        file_size = os.path.getsize(to_keep_file)
        logging.info(f"Keep list file size: {file_size} bytes")
        
        # Check for problematic characters
        with open(to_keep_file, 'rb') as f:
            binary_content = f.read()
            if b'\x00' in binary_content:
                logging.error("ERROR: Keep list file contains null bytes!")
                return

    tree_path = output_file.replace(".fasta", ".tree")

    # Step 2: filter sequences with --seq-k (Stockholm in, Stockholm out)
    if verbose:
        logging.info("Filtering sequences with esl-alimanip --seq-k...")
    tmp_filtered_sto = os.path.join(output_dir, "temp", "filtered.sto")
    
    try:
        subprocess.run(
            [str(alimanip_path), "--seq-k", to_keep_file, "-o", tmp_filtered_sto, tmp_sto],
            check=True
        )
    except subprocess.CalledProcessError as e:
        logging.error(f"esl-alimanip failed: {e}")
        logging.error("This might be due to sequence names mismatch.")
        logging.error(f"First 5 alignment sequence names: {list(full_names)[:5]}")
        logging.error(f"First 5 keep list sequences: {keep[:5]}")
        raise

    # Step 3: reorder to tree order and save Newick tree (incompatible with --seq-k)
    if verbose:
        logging.info("Building tree and reordering sequences...")
    tmp_tree_sto = os.path.join(output_dir, "temp", "tree.sto")
    subprocess.run(
        [str(alimanip_path), "--tree", tree_path, "-o", tmp_tree_sto, tmp_filtered_sto],
        check=True
    )

    # Step 4: convert final Stockholm back to aligned fasta
    if verbose:
        logging.info("Converting final Stockholm back to FASTA format...")
    with open(output_file, "w") as out:
        subprocess.run(
            [str(reformat_path), "afa", tmp_tree_sto],
            stdout=out, check=True
        )

    if verbose:
        logging.info(f"Success! Output written to {output_file}")
        logging.info(f"Debug temp files retained:")
        logging.info(f"  Stockholm input:    {tmp_sto}")
        logging.info(f"  Filtered Stockholm: {tmp_filtered_sto}")
        logging.info(f"  Tree-ordered:       {tmp_tree_sto}")
        logging.info(f"  Keep list:          {to_keep_file}")
    else:
        logging.info("Cleaning up temporary files...")
        shutil.rmtree(os.path.join(output_dir, "temp"))

    # ============================================================================
    # BRANCH COVERAGE ANALYSIS - Calculate per-representative coverage
    # ============================================================================
    logging.info("\n" + "="*70)
    logging.info("BRANCH COVERAGE SUMMARY")
    logging.info("="*70)
    logging.info("(Coverage = kept sequences / (kept + sequences removed by PID to that kept seq))")
    logging.info("="*70)
    
    branch_coverage = {}
    rep_level_coverage = {}  # Per-representative coverage
    
    for branch in sorted(branch_stats.keys()):
        stats = branch_stats[branch]
        kept_seqs = stats['kept']
        removed_seqs = stats['removed']
        
        # For this branch, calculate coverage per kept sequence
        branch_coverage_list = []
        total_with_reps = 0
        total_removed_by_pid = 0
        
        for rep_seq in kept_seqs:
            # Count sequences removed relative to this representative
            removed_by_this_rep = len(pid_pairs.get(rep_seq, []))
            total_removed_by_pid += removed_by_this_rep
            # Total for this rep = itself + sequences removed by PID
            total_for_rep = 1 + removed_by_this_rep
            total_with_reps += total_for_rep
            branch_coverage_list.append({
                'sequence': rep_seq,
                'removed_by_pid': removed_by_this_rep,
                'cluster_size': total_for_rep
            })
        
        # Summary for branch
        num_reps = len(kept_seqs)
        coverage = (num_reps / total_with_reps * 100) if total_with_reps > 0 else 0
        
        branch_coverage[branch] = {
            'num_representatives': num_reps,
            'total_removed_by_pid_to_reps': total_removed_by_pid,
            'total_cluster_size': total_with_reps,
            'coverage_percent': coverage,
            'representative_details': branch_coverage_list
        }
        
        logging.info(f"{branch:15s}: {num_reps:4d}/{total_with_reps:4d} sequences ({coverage:6.2f}% coverage)")
        if verbose:
            logging.debug(f"  Details: {num_reps} reps + {total_removed_by_pid} removed by PID")
            for detail in branch_coverage_list[:3]:  # Show first 3 reps
                logging.debug(f"    {detail['sequence']}: 1 rep + {detail['removed_by_pid']} removed = {detail['cluster_size']} total")
    
    # Write branch coverage to JSON file for visualization
    coverage_file = os.path.join(output_dir, "temp", "branch_coverage.json")
    with open(coverage_file, "w") as f:
        json.dump(branch_coverage, f, indent=2)
    if verbose:
        logging.info(f"\nBranch coverage data written to {coverage_file}")
    
    # Write always-keep reference info
    always_keep_file = os.path.join(output_dir, "temp", "always_keep_references.json")
    with open(always_keep_file, "w") as f:
        json.dump({
            'found': list(always_keep_found),
            'SSL3': always_keep_SSL3,
            'SSL7': always_keep_SSL7,
            'SSL11': always_keep_SSL11
        }, f, indent=2)
    
    if verbose:
        logging.info(f"Always-keep reference data written to {always_keep_file}")
    
    logging.info("="*70 + "\n")
    
    return {
        'branch_coverage': branch_coverage,
        'tree_file': tree_path,
        'output_file': output_file,
        'always_keep_found': list(always_keep_found)
    }



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove sequences from an MSA based on pairwise ID")
    parser.add_argument("--tool-root", required=True, help="Root directory containing EASEL tools")
    parser.add_argument("--maxid", default=0.9, type=float, help="Maximum pairwise ID")
    parser.add_argument("--msafile", required=True, type=str, help="Input MSA file")
    parser.add_argument("-o", "--output-file", required=True, type=str, help="Output MSA file")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    args = parser.parse_args()

    setup_logging(args.verbose, output_dir=os.path.dirname(args.output_file))

    result = main(args.tool_root, args.maxid, args.msafile, args.output_file, verbose=args.verbose)