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
    for line in result.stdout.splitlines():
        if line.startswith("#"): continue
        fields = line.split()
        seq1, seq2, pid = fields[0], fields[1], float(fields[2])
        if pid >= maxid * 100.0:
            if seq1 not in removed:
                removed.add(seq2)   # keep seq1, drop seq2


    # Write keep-list to a temp file and call esl-alimanip --seq-k
    all_seqs_result = subprocess.run(
        [str(alipid_path), "--noheader", msafile],
        capture_output=True, text=True, check=True
    )

    # Collect all sequence names from the alignment
    all_seqs = set()
    for line in result.stdout.splitlines():
        if line.startswith("#"): continue
        f = line.split()
        all_seqs.add(f[0]); all_seqs.add(f[1])

    keep = [s for s in all_seqs if s not in removed]

    with tempfile.NamedTemporaryFile(mode="w", suffix=".keeplist", delete=False) as f:
        tmpfile = f.name
        f.write("\n".join(keep) + "\n")

    output_dir = os.path.dirname(output_file)
    tree_path = os.path.join(output_dir, f"{output_file}.tree")
    alimaniop_path = os.path.join(tool_root, "easel", "miniapps", "esl-alimanip")
    subprocess.run(
        [str(alimaniop_path), "--informat", "afa", "--amino", "--seq-k", 
        tmpfile, msafile, "-o", output_file, "--tree", tree_path],
        check=True
    )
    os.unlink(tmpfile)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove sequences from an MSA based on pairwise ID")
    parser.add_argument("--tool-root", required=True, help="Root directory containing EASEL tools")
    parser.add_argument("--maxid", default=0.9, type=float, help="Maximum pairwise ID")
    parser.add_argument("--msafile", required=True, type=str, help="Input MSA file")
    parser.add_argument("-o","--output-file", required=True, type=str, help="Output MSA file")
    args = parser.parse_args()

    main(args.tool_root, args.maxid, args.msafile, args.output_file)