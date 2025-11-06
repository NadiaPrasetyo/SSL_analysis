#!/usr/bin/env python3
"""
run_jackhmmer_batch.py
-----------------------
Automates running jackhmmer searches for multiple query files against a sequence database.

Usage:
    python run_jackhmmer_batch.py \
        --query-files queries/*.fasta \
        --seq-database uniprot_sprot.fasta \
        --output-dir jackhmmer_results \
        --verbose
"""

import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed


def setup_logging(verbose: bool, output_dir: Path):
    """Configure logging to console and file."""
    log_level = logging.DEBUG if verbose else logging.INFO
    if verbose:
        log_file = output_dir / "jackhmmer_batch.log"

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

def run_jackhmmer(query_file: Path, seq_db: Path, output_dir: Path, verbose: bool = False) -> Path:
    """Run jackhmmer for a single query file."""
    query_name = query_file.stem
    result_file = output_dir / f"{query_name}.jackhmmer.out"
    tblout_file = output_dir / f"{query_name}.jackhmmer.tblout"
    alignment_file = output_dir / f"{query_name}.jackhmmer.sto"

    cmd = [
        "jackhmmer",
        "-N", "5",  # max 5 iterations (default, but explicit)
        "--cpu", "4",  # adjust as needed
        "--tblout", str(tblout_file),
        "-o", str(result_file),
        "--noali",
        "-A", str(alignment_file),
        str(query_file),
        str(seq_db)
    ]

    if verbose:
        logging.debug("Running command: %s", " ".join(cmd))

    try:
        subprocess.run(cmd, check=True, capture_output=not verbose)
        logging.info("[OK] Completed jackhmmer for %s → %s", query_name, result_file)
    except subprocess.CalledProcessError as e:
        logging.error("[FAIL] jackhmmer failed for %s: %s", query_name, e)
        if e.stdout:
            logging.debug("STDOUT:\n%s", e.stdout.decode(errors='ignore'))
        if e.stderr:
            logging.debug("STDERR:\n%s", e.stderr.decode(errors='ignore'))
    except Exception as e:
        logging.exception("Unexpected error running jackhmmer for %s: %s", query_name, e)

    return result_file


def main():
    parser = argparse.ArgumentParser(
        description="Run jackhmmer iteratively on multiple query FASTA files against a sequence database."
    )
    parser.add_argument("--query-files", nargs="+", required=True,
                        help="List of query FASTA files or glob pattern.")
    parser.add_argument("--seq-database", required=True,
                        help="FASTA file containing the target sequence database.")
    parser.add_argument("--output-dir", default="jackhmmer_results",
                        help="Directory to store results (default: jackhmmer_results).")
    parser.add_argument("--verbose", action="store_true",
                        help="Enable verbose progress and debug logging.")

    args = parser.parse_args()

    # Prepare directories
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Setup logging
    setup_logging(args.verbose, output_dir)

    # Verify inputs
    seq_db = Path(args.seq_database)
    if not seq_db.exists():
        logging.error("Sequence database not found: %s", seq_db)
        sys.exit(1)

    query_files = [Path(q) for q in args.query_files]
    query_files = [q for q in query_files if q.exists()]
    if not query_files:
        logging.error("No valid query files found.")
        sys.exit(1)

    logging.info("Found %d query files to process.", len(query_files))
    logging.info("Sequence database: %s", seq_db)
    logging.info("Output directory: %s", output_dir)

    # Run jackhmmer on all queries in parallel
    futures = []
    with ThreadPoolExecutor(max_workers=min(4, len(query_files))) as executor:
        for qfile in query_files:
            futures.append(executor.submit(run_jackhmmer, qfile, seq_db, output_dir, args.verbose))

        for future in as_completed(futures):
            _ = future.result()

    logging.info("All jackhmmer searches completed.")
    logging.info("Results saved in: %s", output_dir)


if __name__ == "__main__":
    main()
