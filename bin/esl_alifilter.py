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

def get_seqs_in_fasta(fasta_file):
    """Extract all sequence names from a FASTA file"""
    seqs = set()
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Remove '>' and take only the first word (in case of descriptions)
                header = line[1:].split()[0]
                seqs.add(header)
    return seqs

def get_seqs_in_stockholm(stockholm_file):
    """Extract all sequence names from a Stockholm format file"""
    seqs = set()
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
                seqs.add(parts[0])
    return seqs
        
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

    # Get all sequences actually in the Stockholm file (canonical names)
    seqs_in_alignment = get_seqs_in_stockholm(tmp_sto)
    if verbose:
        logging.info(f"Found {len(seqs_in_alignment)} sequences in Stockholm alignment")
        logging.debug(f"Sample sequence names: {list(seqs_in_alignment)[:5]}")

    # Run esl-alipid to get all pairwise IDs
    if verbose:
        logging.info("Running esl-alipid to get pairwise identities...")
    result = subprocess.run(
        [str(alipid_path), "--amino", "--informat", "afa", msafile],
        capture_output=True, text=True, check=True
    )

    # Parse pairwise IDs; greedily remove the "worse" of any pair above threshold
    # with preference given to newer sequences
    removed = set()
    all_seqs_ordered = []
    all_seqs_seen = set()
    skipped_count = 0

    for line in result.stdout.splitlines():
        if line.startswith("#"):
            continue
        fields = line.split()
        seq1, seq2, pid = fields[0], fields[1], float(fields[2])
        
        # Only process if both sequences exist in the alignment
        if seq1 not in seqs_in_alignment or seq2 not in seqs_in_alignment:
            if verbose:
                logging.debug(f"Skipping pair: {seq1} vs {seq2} (one or both not in alignment)")
            skipped_count += 1
            continue
        
        for s in (seq1, seq2):
            if s not in all_seqs_seen:
                all_seqs_ordered.append(s)
                all_seqs_seen.add(s)
        
        if pid >= maxid * 100.0:
            year1 = get_year(seq1)
            year2 = get_year(seq2)

            if year1 > 2016 and year2 <= 2016:
                to_remove = seq2  # prefer seq1 (newer)
            elif year2 > 2016 and year1 <= 2016:
                to_remove = seq1  # prefer seq2 (newer)
            else:
                to_remove = seq2  # both same era: fall back to original behaviour

            if to_remove not in removed:
                other = seq1 if to_remove == seq2 else seq2
                if other not in removed:
                    removed.add(to_remove)
                    if verbose:
                        logging.debug(f"Marking for removal: {to_remove} (PID: {pid:.1f}% with {other})")

    if verbose:
        logging.info(f"Skipped {skipped_count} pairs where sequences weren't in alignment")

    keep = [s for s in all_seqs_ordered if s not in removed]

    always_keep_SSL3 = [
        "SSL3|CC1", "SSL3|CC5", "SSL3|CC8", "SSL3|CC22", "SSL3|CC30", "SSL3|CC93"
    ]
    always_keep_SSL7 = [
        "SSL7|CC1", "SSL7|CC5", "SSL7|CC8", "SSL7|CC22", "SSL7|CC30", "SSL7|CC93"
    ]
    always_keep_SSL11 = [
        "SSL11|CC1", "SSL11|CC5", "SSL11|CC8", "SSL11|CC22", "SSL11|CC30", "SSL11|CC93"
    ]

    # Add specific sequences to keep list based on output file name
    # ONLY if they actually exist in the alignment
    added_count = 0
    if "SSL3" in output_file:
        for s in always_keep_SSL3:
            if s in seqs_in_alignment:
                keep.append(s)
                added_count += 1
    elif "SSL7" in output_file:
        for s in always_keep_SSL7:
            if s in seqs_in_alignment:
                keep.append(s)
                added_count += 1
    elif "SSL11" in output_file:
        for s in always_keep_SSL11:
            if s in seqs_in_alignment:
                keep.append(s)
                added_count += 1
    
    if verbose:
        logging.info(f"Added {added_count} reference sequences to keep list")

    keep = list(set(keep))  # deduplicate
    
    if verbose:
        logging.info(f"Final keep list: {len(keep)} sequences")
        if len(keep) == 0:
            logging.warning("WARNING: Keep list is empty! No sequences will be retained.")

    # Write keep-list to a persistent intermediate file (not temp)
    to_keep_file = os.path.join(output_dir, "temp", "to_keep.txt")
    with open(to_keep_file, "w") as f:
        f.write("\n".join(keep) + "\n")
    
    if verbose:
        logging.info(f"Wrote keep list to {to_keep_file}")

    tree_path = output_file.replace(".fasta", ".tree")

    # Step 2: filter sequences with --seq-k (Stockholm in, Stockholm out)
    if verbose:
        logging.info("Filtering sequences with esl-alimanip --seq-k...")
    tmp_filtered_sto = os.path.join(output_dir, "temp", "filtered.sto")
    subprocess.run(
        [str(alimanip_path), "--seq-k", to_keep_file, "-o", tmp_filtered_sto, tmp_sto],
        check=True
    )

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



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove sequences from an MSA based on pairwise ID")
    parser.add_argument("--tool-root", required=True, help="Root directory containing EASEL tools")
    parser.add_argument("--maxid", default=0.9, type=float, help="Maximum pairwise ID")
    parser.add_argument("--msafile", required=True, type=str, help="Input MSA file")
    parser.add_argument("-o", "--output-file", required=True, type=str, help="Output MSA file")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    args = parser.parse_args()

    setup_logging(args.verbose, output_dir=os.path.dirname(args.output_file))

    main(args.tool_root, args.maxid, args.msafile, args.output_file, verbose=args.verbose)