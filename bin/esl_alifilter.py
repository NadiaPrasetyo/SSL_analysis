#!/usr/bin/env python3
import sys
import subprocess
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



def get_seqs_in_stockholm(stockholm_file):
    """Extract all sequence names from a Stockholm format file.
    Returns both full names (with prefix) and base names (without prefix)."""
    full_names = set()
    base_names = defaultdict(list)
    
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

                base_names[base_name].append(full_name)

    return full_names, base_names

def cluster_sequences_by_pid(pairs, full_names, base_names_map, maxid, always_keep_sequences, verbose=False):
    representative = {}
    represented = defaultdict(set)

    # Initialize each sequence as its own cluster
    for seq in full_names:
        representative[seq] = seq
        represented[seq].add(seq)

    def score(seq):
        if seq in always_keep_sequences:
            return (2, 0)
        return (1, -seq.count("unknown"))

    # --- CLUSTERING ---
    for seq1, seq2, pid in pairs:
        if pid < maxid * 100.0:
            break

        full_seqs1 = base_names_map.get(seq1, [])
        full_seqs2 = base_names_map.get(seq2, [])

        if not full_seqs1 or not full_seqs2:
            continue

        # 🔥 THIS is where your snippet goes
        for f1 in full_seqs1:
            for f2 in full_seqs2:

                rep1 = representative[f1]
                rep2 = representative[f2]

                if rep1 == rep2:
                    continue

                if score(rep1) >= score(rep2):
                    keep_rep, remove_rep = rep1, rep2
                else:
                    keep_rep, remove_rep = rep2, rep1

                # Merge clusters safely
                if remove_rep not in represented:
                    continue

                for s in represented[remove_rep]:
                    representative[s] = keep_rep
                    represented[keep_rep].add(s)

                del represented[remove_rep]

    # --- Representatives ---
    keep = sorted(represented.keys())

    # --- Enforce always-keep sequences ---
    for base in always_keep_sequences:
        full_list = base_names_map.get(base, [])

        for full_s in full_list:
            if full_s not in representative:
                continue

            rep = representative[full_s]

            if rep != full_s:
                if rep in represented and full_s in represented[rep]:
                    represented[rep].remove(full_s)

                represented[full_s] = {full_s}
                representative[full_s] = full_s
                keep.append(full_s)

    keep = sorted(set(keep))

    # --- Coverage ---
    total_original_sequences = len(full_names)

    coverage_report = {}
    for rep, cluster in represented.items():
        num = len(cluster)
        coverage_percent = (num / total_original_sequences) * 100.0

        coverage_report[rep] = {
            "num": num,
            "coverage_percent": round(coverage_percent, 2),
            "removed_sequences": sorted([s for s in cluster if s != rep])
        }

    total_coverage = sum(v["num"] for v in coverage_report.values())
    assert total_coverage == total_original_sequences

    if verbose:
        logging.info(f"Clusters formed: {len(keep)}")

    return keep, coverage_report


def main(tool_root, maxid, msafile, output_file, verbose=False):
    alipid_path = os.path.join(tool_root, "easel", "miniapps", "esl-alipid")
    alimanip_path = os.path.join(tool_root, "easel", "miniapps", "esl-alimanip")
    reformat_path = os.path.join(tool_root, "easel", "miniapps", "esl-reformat")

    # Step 1: convert input AFA -> Stockholm first (to get canonical sequence names)
    output_dir = os.path.dirname(output_file) or "."
    os.makedirs(os.path.join(output_dir, "temp"), exist_ok=True)

    tmp_sto = os.path.join(output_dir, "temp", "tmp.sto")

    if verbose:
        logging.info(f"Converting {msafile} to Stockholm...")

    with open(tmp_sto, "w") as out:
        subprocess.run(
            [reformat_path, "--informat", "afa", "stockholm", msafile],
            stdout=out, check=True
        )

    # Get all sequences actually in the Stockholm file
    full_names, base_names_map = get_seqs_in_stockholm(tmp_sto)

    if verbose:
        logging.info(f"Sequences in alignment: {len(full_names)}")

    if verbose:
        logging.info("Running esl-alipid...")

    result = subprocess.run(
        [alipid_path, "--amino", "--informat", "afa", msafile],
        capture_output=True, text=True, check=True
    )

    # --- Parse + clean PID pairs ---
    pairs = []
    pair_best_pid = {}

    for line in result.stdout.splitlines():
        if line.startswith("#"):
            continue

        fields = line.split()
        if len(fields) < 3:
            continue

        seq1, seq2, pid = fields[0], fields[1], float(fields[2])

        if seq1 == seq2: #skip self
            continue

        key = tuple(sorted((seq1, seq2)))
        pair_best_pid[key] = max(pid, pair_best_pid.get(key, 0.0))

    pairs = [(k[0], k[1], v) for k, v in pair_best_pid.items()]
    pairs.sort(key=lambda x: x[2], reverse=True)

    if verbose:
        logging.info(f"Unique PID pairs: {len(pairs)}")

    # --- Always keep sets ---
    always_keep = set()
    if "SSL3" in output_file:
        always_keep = {"SSL3|CC1","SSL3|CC5","SSL3|CC8","SSL3|CC22","SSL3|CC30","SSL3|CC93"}
    elif "SSL7" in output_file:
        always_keep = {"SSL7|CC1","SSL7|CC5","SSL7|CC8","SSL7|CC22","SSL7|CC30","SSL7|CC93"}
    elif "SSL11" in output_file:
        always_keep = {"SSL11|CC1","SSL11|CC5","SSL11|CC8","SSL11|CC22","SSL11|CC30","SSL11|CC93"}

    # --- Cluster ---
    keep, coverage_report = cluster_sequences_by_pid(
        pairs,
        full_names,
        base_names_map,
        maxid,
        always_keep,
        verbose=verbose
    )

    if verbose:
        logging.info(f"Final keep count: {len(keep)}")

    # --- Write keep list ---
    to_keep_file = os.path.join(output_dir, "temp", "to_keep.txt")
    with open(to_keep_file, "w") as f:
        for seq in keep:
            f.write(seq + "\n")

    # --- Coverage output ---
    if verbose:
        cov_file = os.path.join(output_dir, "temp", "sequence_coverage.json")
        with open(cov_file, "w") as f:
            json.dump(coverage_report, f, indent=2)

        logging.info(f"Coverage written to {cov_file}")

    # --- Filter alignment ---
    tmp_filtered = os.path.join(output_dir, "temp", "filtered.sto")
    subprocess.run(
        [alimanip_path, "--seq-k", to_keep_file, "-o", tmp_filtered, tmp_sto],
        check=True
    )

    # --- Build tree ---
    tree_path = output_file.replace(".fasta", ".tree")
    tmp_tree = os.path.join(output_dir, "temp", "tree.sto")

    subprocess.run(
        [alimanip_path, "--tree", tree_path, "-o", tmp_tree, tmp_filtered],
        check=True
    )

    # --- Convert back to FASTA ---
    with open(output_file, "w") as out:
        subprocess.run(
            [reformat_path, "afa", tmp_tree],
            stdout=out, check=True
        )

    if not verbose:
        shutil.rmtree(os.path.join(output_dir, "temp"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove sequences from an MSA based on pairwise ID")
    parser.add_argument("--tool-root", required=True, help="Root directory containing EASEL tools")
    parser.add_argument("--maxid", default=0.9, type=float, help="Maximum pairwise ID")
    parser.add_argument("--msafile", required=True, type=str, help="Input MSA file")
    parser.add_argument("-o", "--output-file", required=True, type=str, help="Output MSA file")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    args = parser.parse_args()

    setup_logging(args.verbose, os.path.dirname(args.output_file))
    main(args.tool_root, args.maxid, args.msafile, args.output_file, args.verbose)