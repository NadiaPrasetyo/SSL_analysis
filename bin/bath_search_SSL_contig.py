import sys
import os
import argparse
import logging
from pathlib import Path
import tempfile
import subprocess

def setup_logging(verbose: bool, output_dir: Path):
    """Configure logging to console and file."""
    log_level = logging.DEBUG if verbose else logging.INFO
    log_file = output_dir / "fetch_pubMLST_contigs.log"
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file, mode='a') if verbose else logging.NullHandler(),
            logging.StreamHandler(sys.stdout)
        ]
    )
    if verbose:
        logging.info("Logging initialized. Log file: %s", log_file)


def check_tool(tool_path):
    if not os.path.exists(tool_path):
        logging.error(f"Tool not found at {tool_path}")
        sys.exit(1)
    if not os.access(tool_path, os.X_OK):
        logging.error(f"Tool at {tool_path} is not executable")
        sys.exit(1)


def already_done(output_dir: Path, file_stem: str) -> bool:
    """Return True if both expected output files exist and are non-empty."""
    out_file = output_dir / f"{file_stem}_bath_results.out"
    tbl_file = output_dir / f"{file_stem}_bath_results.tbl"
    return out_file.exists() and out_file.stat().st_size > 0 \
       and tbl_file.exists() and tbl_file.stat().st_size > 0


def main():
    parser = argparse.ArgumentParser(description="Batch search SSL contig")
    parser.add_argument("-i", "--input", type=str, required=True, help="Input Query FASTA")
    parser.add_argument("-t", "--target", type=str, required=True, help="Target Directory that contains FASTA files")
    parser.add_argument("-o", "--output", type=str, default="bath_results", help="Output folder (default: bath_results)")
    parser.add_argument("--tool-path", type=str, required=True, help="Path to the tool root directory")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    parser.add_argument("--reset", action="store_true", help="Reprocess all files, ignoring existing outputs")
    args = parser.parse_args()

    # Set up logging
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(args.verbose, output_dir)

    tool_path = args.tool_path
    bathsearch_path = os.path.join(tool_path, "bathsearch")
    bathbuild_path = os.path.join(tool_path, "bathbuild")
    check_tool(bathbuild_path)
    check_tool(bathsearch_path)

    target_files = sorted([f for f in os.listdir(args.target) if f.endswith(".fasta")])
    if not target_files:
        logging.error("No FASTA files found in the target directory.")
        sys.exit(1)

    # ── bathbuild ──────────────────────────────────────────────────────────────
    temp_dir = tempfile.mkdtemp()
    # Use only the stem of the input filename to avoid nested subdirs in temp
    input_stem = Path(args.input).stem
    bhmm_file = os.path.join(temp_dir, f"{input_stem}.bhmm")

    logging.info("Running bathbuild to create HMM profile…")
    cmd = [str(bathbuild_path), "--unali", bhmm_file, args.input]
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logging.error("Error running bathbuild command: %s", e)
        sys.exit(1)

    # ── bathsearch loop ────────────────────────────────────────────────────────
    total     = len(target_files)
    skipped   = 0
    processed = 0

    for idx, file in enumerate(target_files, start=1):
        file_stem   = Path(file).stem
        output_file = output_dir / f"{file_stem}_bath_results"

        if not args.reset and already_done(output_dir, file_stem):
            logging.info("[%d/%d] Skipping %s (output already exists).", idx, total, file)
            skipped += 1
            continue

        logging.info("[%d/%d] Processing: %s", idx, total, file)
        cmd = [
            str(bathsearch_path),
            "-o",       f"{output_file}.out",
            "--tblout", f"{output_file}.tbl",
            bhmm_file,
            str(os.path.join(args.target, file))
        ]

        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            logging.error("Error running bathsearch for %s: %s", file, e)
            logging.info("Re-run the same command to continue from this point.")
            sys.exit(1)

        processed += 1
        logging.info("Completed %s (%d/%d).", file, idx, total)

    logging.info("Done. Processed: %d, Skipped: %d, Total: %d.", processed, skipped, total)


if __name__ == "__main__":
    main()