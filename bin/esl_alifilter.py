#!/usr/bin/env python3
# Usage: esl_alifilter.py <maxid> <msafile>
# e.g.:  esl_alifilter.py 0.90 my_alignment.sto
import sys
import subprocess
import tempfile
import os
import argparse

def main(tool_root, maxid, msafile, output_file):
    alipid_path = os.path.join(tool_root, "easel", "miniapps", "esl-alipid")

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

        # Collect all sequence names (BUG FIX 1: was reading from wrong variable
        # all_seqs_result.stdout which was never populated correctly)
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

    # add specific sequences to keep list
    if "SSL3" in output_file:
        keep.extend(always_keep_SSL3)
    elif "SSL7" in output_file:
        keep.extend(always_keep_SSL7)
    elif "SSL11" in output_file:
        keep.extend(always_keep_SSL11)

    keep = list(set(keep))  # deduplicate
    
    # write to keep list on an intermediate file (NOT TEMP)
    to_keep_file = output_file + ".keep"
    with open(to_keep_file, "w") as f:
        f.write("\n".join(keep))


    # the suffix directly — don't re-join with output_dir or the path doubles up
    tree_path = f"{output_file}.tree"

    alimanip_path = os.path.join(tool_root, "easel", "miniapps", "esl-alimanip")  # BUG FIX 3: typo alimaniop -> alimanip

    subprocess.run(
        [str(alimanip_path), "--informat", "afa", "--amino", "--seq-k",
         to_keep_file, "-o", output_file, "--tree", tree_path, msafile],
        check=True
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove sequences from an MSA based on pairwise ID")
    parser.add_argument("--tool-root", required=True, help="Root directory containing EASEL tools")
    parser.add_argument("--maxid", default=0.9, type=float, help="Maximum pairwise ID")
    parser.add_argument("--msafile", required=True, type=str, help="Input MSA file")
    parser.add_argument("-o", "--output-file", required=True, type=str, help="Output MSA file")
    args = parser.parse_args()
    main(args.tool_root, args.maxid, args.msafile, args.output_file)