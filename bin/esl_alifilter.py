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

def get_seqs_in_alignment(msafile):
    """Extract all sequence names from the alignment file"""
    seqs = set()
    with open(msafile, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                seqs.add(line[1:])  # remove '>' prefix
    return seqs
        
def main(tool_root, maxid, msafile, output_file, verbose=False):
    alipid_path      = os.path.join(tool_root, "easel", "miniapps", "esl-alipid")
    alimanip_path    = os.path.join(tool_root, "easel", "miniapps", "esl-alimanip")
    reformat_path = os.path.join(tool_root, "easel", "miniapps", "esl-reformat")
    output_dir = os.path.dirname(output_file)

    # Get all sequences actually in the alignment
    seqs_in_alignment = get_seqs_in_alignment(msafile)
    if verbose:
        logging.info(f"Found {len(seqs_in_alignment)} sequences in alignment")

    # Run esl-alipid to get all pairwise IDs
    result = subprocess.run(
        [str(alipid_path), "--amino", "--informat", "afa", msafile],
        capture_output=True, text=True, check=True
    )

    # Parse pairwise IDs; greedily remove the "worse" of any pair above threshold
    # with preference given to newer sequences
    removed = set()
    all_seqs_ordered = []
    all_seqs_seen = set()

    for line in result.stdout.splitlines():
        if line.startswith("#"):
            continue
        fields = line.split()
        seq1, seq2, pid = fields[0], fields[1], float(fields[2])
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
    if "SSL3" in output_file:
        keep.extend([s for s in always_keep_SSL3 if s in seqs_in_alignment])
    elif "SSL7" in output_file:
        keep.extend([s for s in always_keep_SSL7 if s in seqs_in_alignment])
    elif "SSL11" in output_file:
        keep.extend([s for s in always_keep_SSL11 if s in seqs_in_alignment])

    keep = list(set(keep))  # deduplicate

    # Filter keep list to only include sequences that actually exist
    keep = [s for s in keep if s in seqs_in_alignment]
    
    if verbose:
        logging.info(f"After filtering: keeping {len(keep)} sequences")
        invalid = [s for s in list(set(keep)) if s not in seqs_in_alignment]
        if invalid:
            logging.warning(f"Removed {len(invalid)} invalid sequences from keep list")

    # Write keep-list to a persistent intermediate file (not temp)
    os.makedirs(os.path.join(output_dir, "temp"), exist_ok=True)
    to_keep_file = os.path.join(output_dir, "temp", "to_keep.txt")
    with open(to_keep_file, "w") as f:
        f.write("\n".join(keep) + "\n")

    tree_path = output_file.replace(".fasta", ".tree")

    # Step 1: convert input AFA -> Stockholm (--seq-k requires Stockholm format)
    tmp_sto = os.path.join(output_dir, "temp", "tmp.sto")
    with open(tmp_sto, "w") as out:
        subprocess.run(
            [str(reformat_path), "--informat", "afa", "stockholm", msafile],
            stdout=out, check=True
        )

    # Step 2: filter sequences with --seq-k (Stockholm in, Stockholm out)
    tmp_filtered_sto = os.path.join(output_dir, "temp", "filtered.sto")
    subprocess.run(
        [str(alimanip_path), "--seq-k", to_keep_file, "-o", tmp_filtered_sto, tmp_sto],
        check=True
    )

    # Step 3: reorder to tree order and save Newick tree (incompatible with --seq-k)
    tmp_tree_sto = os.path.join(output_dir, "temp", "tree.sto")
    subprocess.run(
        [str(alimanip_path), "--tree", tree_path, "-o", tmp_tree_sto, tmp_filtered_sto],
        check=True
    )

    # Step 4: convert final Stockholm back to aligned fasta
    with open(output_file, "w") as out:
        subprocess.run(
            [str(reformat_path), "afa", tmp_tree_sto],
            stdout=out, check=True
        )

    if verbose:
        logging.info(f"Debug temp files retained:")
        logging.info(f"  Stockholm input:    {tmp_sto}")
        logging.info(f"  Filtered Stockholm: {tmp_filtered_sto}")
        logging.info(f"  Tree-ordered:       {tmp_tree_sto}")
        logging.info(f"  Keep list:          {to_keep_file}")
    else:
        logging.info("Debug temp files removed.")
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