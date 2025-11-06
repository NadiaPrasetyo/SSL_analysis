#!/usr/bin/env python3
"""
Script: compile_6frame_proteins.py
Description: Concatenate all *_6frame_proteins.fasta files in specified directories
             into one combined FASTA file (genomes_6frame_proteins.fasta).

Usage:
    python compile_6frame_proteins.py \
        --base-dir data/my_genomes \
        --input-dirs genome1 genome2 genome3 \
        --output-dir results \
        --verbose
"""

import os
import sys
import argparse
import logging
from pathlib import Path

def setup_logging(output_dir: Path, verbose: bool):
    """Set up logging to file and console."""
    log_format = "%(asctime)s [%(levelname)s] %(message)s"
    if verbose:
        log_file = output_dir / "compile_6frame_proteins.log"
    logging.basicConfig(
        level= logging.INFO,
        format=log_format,
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_file, mode='w') if verbose else logging.NullHandler()
        ]
    )

def collect_fasta_files(base_dir: Path, input_dirs: list):
    """Collect all *_6frame_proteins.fasta files from the input directories."""
    fasta_files = []
    for d in input_dirs:
        dir_path = base_dir / d
        if not dir_path.exists():
            logging.warning(f"Directory not found: {dir_path}")
            continue
        matches = list(dir_path.glob("*_6frame_proteins.fasta"))
        if not matches:
            logging.warning(f"No matching FASTA files in {dir_path}")
        else:
            fasta_files.extend(matches)
            for f in matches:
                logging.debug(f"Found: {f}")
    return fasta_files

def concatenate_fastas(fasta_files: list, output_path: Path):
    """Concatenate all FASTA files into one."""
    with open(output_path, 'w') as outfile:
        for f in fasta_files:
            logging.info(f"Adding: {f}")
            with open(f, 'r') as infile:
                for line in infile:
                    outfile.write(line)
    logging.info(f"Combined FASTA written to: {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Compile all *_6frame_proteins.fasta files.")
    parser.add_argument("--base-dir", required=True, help="Base directory under data/ containing genome dirs")
    parser.add_argument("--input-dirs", required=True, nargs='+', help="List of directories containing *_6frame_proteins.fasta")
    parser.add_argument("--output-dir", help="Directory to save the final FASTA (default: base-dir)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")

    args = parser.parse_args()

    base_dir = Path("data") / args.base_dir
    output_dir = Path(args.output_dir) if args.output_dir else base_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    

    setup_logging(output_dir, args.verbose)
    logging.info(f"Base directory: {base_dir}")
    logging.info(f"Input directories: {', '.join(args.input_dirs)}")
    logging.info(f"Output directory: {output_dir}")

    fasta_files = collect_fasta_files(base_dir, args.input_dirs)
    if not fasta_files:
        logging.error("No *_6frame_proteins.fasta files found. Exiting.")
        sys.exit(1)

    output_fasta = output_dir / "genomes_6frame_proteins.fasta"
    concatenate_fastas(fasta_files, output_fasta)

if __name__ == "__main__":
    main()
