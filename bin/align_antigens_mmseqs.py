"""
align_antigens_mmseqs.py

Command-line tool to align antigen sequences against strain genome sequences using MMseqs2.

Overview:
    - Converts a compiled antigen CSV file to FASTA format.
    - Aligns antigens to each strain's genome or proteome using MMseqs2.
    - Extracts the best alignment hits for each antigen.
    - Saves results as TSV (all alignments, best hits) and FASTA (matched antigen regions).

Arguments:
    pathogen_directory (str): Subdirectory under `data/` containing pathogen data.
    pathogen_name (str): Scientific name of the pathogen (used to infer antigen file names).
    --genome-dir (str, optional): Subdirectory under pathogen_directory containing strain genome or proteome FASTA files (default: strain_genomes).
    --threads (int, optional): Number of parallel workers (default: 4).
    --mode (str, optional): Alignment mode, either "protein" or "nucleotide" (default: protein).
    --output-dir (str, optional): Output subdirectory (default: mmseqs_results).
    --fetch-qseq (flag, optional): If set, includes query sequences in MMseqs2 output.
    --verbose (flag, optional): Enables verbose logging to console and file.

Requirements:
    - MMseqs2 installed and available in PATH.
    - Antigen CSV file with columns: uniprot_accession, protein_name, sequence (and nucleotide_sequence if in nucleotide mode).
    - Strain genome or proteome FASTA files named appropriately.

Usage Example:
    python align_antigens_mmseqs.py sars_cov_2 "SARS-CoV-2" --threads 8 --fetch-qseq --mode protein

Outputs:
    <output_dir>/<strain>_alignment.tsv         # Raw MMseqs2 alignments
    <output_dir>/<strain>_best_hits.tsv         # Best hit per antigen
    <output_dir>/<strain>_matched_antigens.fasta # FASTA of matched antigen regions

Author: Nadia
"""

#!/usr/bin/env python3
import os
import csv
import sys
import subprocess
import tempfile
import shutil
import argparse
from Bio import SeqIO
from pathlib import Path
from pathlib import Path
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed

def setup_logging(verbose, output_dir):
    """
    Configures logging to output to both console and a log file.
    Args:
        verbose (bool): If True, enables verbose logging.
        output_dir (str): Path to the output directory.
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    log_file = Path(output_dir) / "mmseqs_alignment.log"
    os.makedirs(output_dir, exist_ok=True)
    
    handlers = [logging.StreamHandler()]  # Log to console
    if verbose:
        handlers.append(logging.FileHandler(log_file))  # Log to file

    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=handlers
    )

def extract_antigens_to_fasta(csv_path, fasta_path, mode="protein"):
    """
    Converts antigen CSV records into a FASTA file.
    Each row in the input CSV is formatted into a FASTA record using the UniProt accession,
    protein name, and sequence fields.
    Parameters:
        csv_path (str): Path to the antigen CSV file.
        fasta_path (str): Output path for the generated FASTA file.
        mode (str): Mode for extraction ("protein" or "nucleotide").
    Returns:
        None
    """
    with open(csv_path, newline='') as csvfile, open(fasta_path, 'w') as f_out:
        reader = csv.DictReader(csvfile)
        for idx, row in enumerate(reader):
            acc = row['uniprot_accession']
            name = row['protein_name']
            if mode == "protein":
                seq = row['sequence'].replace('\r', '').replace('\n', '')
            elif mode == "nucleotide":
                seq = row['nucleotide_sequence'].replace('\r', '').replace('\n', '')
            f_out.write(f">antigen_{idx}|{acc}|{name}\n{seq}\n")

def run_mmseqs2_and_process(strain_fasta_path, antigen_fasta, results_dir, fetch_qseq=False, mode="protein", dir_prefix=""):
    """
    Runs MMseqs2 alignment of antigens against a strain's genome or proteome,
    processes the results to extract best hits, and saves outputs.
    Args:
        strain_fasta_path (str): Path to the strain's genome or proteome FASTA
        antigen_fasta (str): Path to the antigen FASTA file.
        results_dir (Path): Directory to save output files.
        fetch_qseq (bool): Whether to include query sequences in the output.
        mode (str): Alignment mode ("protein" or "nucleotide").
        dir_prefix (str): Optional prefix for output filenames.
    Returns:
        None
    """
    strain_fasta = Path(strain_fasta_path)
    strain_name = strain_fasta.stem.replace("_6frame_proteins", "")
    
    # Add pathogen directory name prefix to outputs
    prefix = f"{dir_prefix}_" if dir_prefix else ""
    raw_result = results_dir / f"{prefix}{strain_name}_alignment.tsv"
    best_result = results_dir / f"{prefix}{strain_name}_best_hits.tsv"
    antigen_seqs_out = results_dir / f"{prefix}{strain_name}_matched_antigens.fasta"
    ...
    extract_best_hits_with_sequences(strain_fasta_path, raw_result, best_result, antigen_seqs_out, fetch_qseq)

    # Remove empty output files
    for path in [raw_result, best_result, antigen_seqs_out]:
        if not path.exists() or os.path.getsize(path) == 0:
            logging.info(f"[⚠] Removing empty file: {path}")
            path.unlink(missing_ok=True)


def extract_best_hits_with_sequences(strain_fasta_path, raw_tsv_path, output_tsv_path, fasta_out_path, fetch_qseq):
    """
    Extracts the best alignment hits from MMseqs2 results and writes them to TSV and FASTA files.
    For each query, retains the best hit based on percent identity and outputs:
    - A TSV summary of the best hits.
    - A FASTA file containing the matched target sequence slices.
    Parameters:
        strain_fasta_path (str): Path to the FASTA file containing all target sequences.
        raw_tsv_path (str): Path to the raw MMseqs2 alignment TSV file.
        output_tsv_path (str): Output path for the filtered best hits TSV file.
        fasta_out_path (str): Output path for the matched target sequences FASTA file.
        fetch_qseq (bool): Whether to include the query sequence in the output TSV.
    Returns:
        None
    """
    best_hits = {}

    # Load all target sequences into a dict for slicing
    target_seqs = {}
    for record in SeqIO.parse(str(strain_fasta_path), "fasta"):
        target_seqs[record.id] = str(record.seq)

    with open(raw_tsv_path) as f:
        lines = f.readlines()
        if lines and lines[0].lower().startswith("query"):
            lines = lines[1:]

        for line in lines:
            parts = line.strip().split('\t')
            if len(parts) < 13:
                continue

            try:
                query, target = parts[0], parts[1]
                pident = float(parts[2])
                evalue = parts[5]
                mismatch = parts[7]
                qcov = parts[8]
                tcov = parts[9]
                tstart, tend = int(parts[10]), int(parts[11])
                taln = parts[12]
                qseq = parts[13] if fetch_qseq and len(parts) > 13 else ""
            except (IndexError, ValueError):
                continue

            if query not in best_hits or pident > best_hits[query]['pident']:
                target_seq = target_seqs.get(target, "")

                # MMseqs uses 1-based inclusive coordinates; Python slicing is 0-based and exclusive at end
                if target_seq:
                    if tstart < tend:
                        tseq_slice = target_seq[tstart - 1:tend]
                    elif tstart > tend:
                        tseq_slice = target_seq[tend - 1:tstart][::-1]  # Reverse if alignment is on opposite strand
                    else:
                        tseq_slice = ""
                else:
                    tseq_slice = ""

                best_hits[query] = {
                    'query': query,
                    'target': target,
                    'pident': pident,
                    'evalue': evalue,
                    'mismatch': mismatch,
                    'qcov': qcov,
                    'tcov': tcov,
                    'tstart': tstart,
                    'tend': tend,
                    'taln': taln,
                    'tseq_slice': tseq_slice,
                    'qseq': qseq,
                }

    headers = ["query", "target", "pident", "evalue", "mismatch", "qcov", "tcov", "tstart", "tend", "taln"]
    if fetch_qseq:
        headers.append("qseq")

    with open(output_tsv_path, 'w') as f_out:
        f_out.write('\t'.join(headers) + '\n')
        for hit in best_hits.values():
            row = [str(hit[h]) for h in headers]
            f_out.write('\t'.join(row) + '\n')

    with open(fasta_out_path, 'w') as fasta_out:
        for hit in best_hits.values():
            header = f"{hit['query']}|{hit['target']}|tpos:{hit['tstart']}-{hit['tend']}"
            fasta_out.write(f">{header}\n{hit['tseq_slice']}\n")

def main(pathogen_dir, pathogen_name, genome_dir, num_threads, output_dir, fetch_qseq, mode):
    """
    Entry point to execute the antigen-to-strain alignment workflow.
    Validates input files and directories, prepares an antigen FASTA file, and runs
    MMseqs2 alignments in parallel for each strain genome.
    Args:
        pathogen_dir (str): Directory name under `data/` containing pathogen data.
        pathogen_name (str): Scientific name of the organism.
        num_threads (int): Number of parallel processes to use for alignment.
        output_dir (str): Name of the output subdirectory for results.
        fetch_qseq (bool): If True, include query sequences in the results.
    Returns:
        None
    """
    base_dir = Path(f"data/{pathogen_dir}")
    pathogen_tag = pathogen_name.replace(" ", "_").lower()
    antigen_csv = base_dir / f"{pathogen_tag}_compiled_proteins.csv"
    strain_dir = base_dir / genome_dir
    results_dir = Path(base_dir/output_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    if mode == "protein":
        strain_files = list(strain_dir.glob("*_6frame_proteins.fasta"))
    else:  # nucleotide
        strain_files = [f for f in strain_dir.glob("*.fasta") if not f.name.endswith("_proteins.fasta")]
        logging.info(f"Running in nucleotide mode. Found {len(strain_files)} files: {strain_files} .")

    if not antigen_csv.exists():
        logging.error(f"Antigen CSV file {antigen_csv} does not exist.")
        sys.exit(1)
    if not strain_dir.exists() or not strain_files:
        logging.error(f"No strain FASTA files found in {strain_dir}.")
        sys.exit(1)
    if not shutil.which("mmseqs"):
        logging.error("MMseqs2 is not installed or not found in PATH.")
        sys.exit(1)
    if not base_dir.exists() or not base_dir.is_dir():
        logging.error(f"Error: Base directory {base_dir} does not exist.")
        sys.exit(1)

    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as tmp:
        antigen_fasta = tmp.name
    extract_antigens_to_fasta(antigen_csv, antigen_fasta, mode)

    logging.info(f"Antigen FASTA file created with mode {mode}: ({antigen_fasta}).")

    strain_fastas = list(strain_files)
    logging.info(f"Running MMseqs2 on {len(strain_fastas)} strains with {num_threads} workers...")

    with ProcessPoolExecutor(max_workers=num_threads) as executor:
        futures = {
            executor.submit(run_mmseqs2_and_process, f, antigen_fasta, results_dir, fetch_qseq, mode, genome_dir): f
            for f in strain_fastas
        }
        for future in as_completed(futures):
            future.result()

    os.remove(antigen_fasta)
    logging.info("All alignments complete.")

if __name__ == "__main__":
    """
    Command-line interface for aligning antigens to strain genomes using MMseqs2.
    Parse command-line arguments and run the main alignment function.
    """
    
    parser = argparse.ArgumentParser(
        description="Align antigens to strain genomes using MMseqs2"
    )
    parser.add_argument("pathogen_directory", help="Directory name under data/")
    parser.add_argument("pathogen_name", help='Prefix used in filenames (e.g., "staphylococcus aureus")')
    parser.add_argument("--genome-dir", help="Subdirectory under pathogen_directory containing strain genome or proteome FASTA files (default: strain_genomes)", default="strain_genomes")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads (default: 4)")
    parser.add_argument("--mode", choices=["protein", "nucleotide"], default="protein", help="Alignment mode (default: protein)")
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Custom output directory (default: data/<pathogen_directory>/mmseqs_results)"
    )
    parser.add_argument(
        "--fetch-qseq",
        action="store_true",
        help="Include qseq (query sequence) in MMseqs2 output"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging to console and file"
    )

    args = parser.parse_args()
    setup_logging(args.verbose, output_dir=args.output_dir if args.output_dir else f"data/{args.pathogen_directory}/mmseqs_results")

    if args.threads < 2:
        logging.error("Please specify at least 2 threads with --threads.")
        sys.exit(1)

    # Set default output dir if not provided
    if args.output_dir is None:
        args.output_dir = f"mmseqs_results"

    main(args.pathogen_directory, args.pathogen_name, args.genome_dir, args.threads, args.output_dir, args.fetch_qseq, args.mode)
