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

def main(tool_root, maxid, msafile, output_file, verbose=False):
    alipid_path      = os.path.join(tool_root, "easel", "miniapps", "esl-alipid")
    alimanip_path    = os.path.join(tool_root, "easel", "miniapps", "esl-alimanip")
    reformat_path = os.path.join(tool_root, "easel", "miniapps", "esl-reformat")
    output_dir = os.path.dirname(output_file)

    # Run esl-alipid to get all pairwise IDs
    result = subprocess.run(
        [str(alipid_path), "--amino", "--informat", "afa", msafile],
        capture_output=True, text=True, check=True
    )

    # Parse pairwise IDs; greedily remove the "worse" of any pair above threshold
    # (keeps first-seen sequence, discards its high-id partners -- simple greedy)
    removed = set()
    all_seqs_ordered = []  # preserve insertion order for keep-list
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
            if seq1 not in removed:
                removed.add(seq2)   # keep seq1, drop seq2

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
    if "SSL3" in output_file:
        keep.extend(always_keep_SSL3)
    elif "SSL7" in output_file:
        keep.extend(always_keep_SSL7)
    elif "SSL11" in output_file:
        keep.extend(always_keep_SSL11)

    keep = list(set(keep))  # deduplicate

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
