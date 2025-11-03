#!/usr/bin/env python3
"""
fetch_ncbi_strain_genomes.py

Rewritten from fetch_NCBI_strain_genome.sh (previously used NCBI Entrez Direct).
Uses NCBI Datasets v2 API for fetching complete genomes.

Dependencies:
    pip install requests tqdm

Features:
    - Fetch random complete genomes for a given taxon.
    - Fetch genomes based on a CSV file containing RefSeq IDs.
    - Automatically downloads and extracts genome FASTA and protein FASTA files.
    - Handles pagination for large datasets.

Usage:
    Fetch random genomes:
        python fetch_ncbi_strain_genomes.py --random "Staphylococcus aureus" \
            --random-num 5 data/staph_aureus

    Fetch genomes from a CSV file:
        python fetch_ncbi_strain_genomes.py data/staph_aureus strains.csv

Arguments:
    pathogen_dir: Output directory under the "data/" folder.
    --csv_file: Optional CSV file with strain names and RefSeq IDs.
    --threads: Number of threads to use (default: 4).
    --random: Fetch random genomes for the specified pathogen.
    --random-num: Number of random genomes to fetch (default: 5).
    --output-dir: Custom output directory (default: pathogen_dir/strain_genomes).
"""

import argparse
import csv
import logging
import os
import random
import shutil
import tempfile
import zipfile
import time
import requests
from tqdm import tqdm
from functools import wraps

# ---------------------------
# Configuration
# ---------------------------
API_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome"
DEFAULT_THREADS = 4
DEFAULT_RANDOM_NUM = 5
RATE_LIMIT_RPS = 5          # 5 requests per second
MIN_REQUEST_INTERVAL = 1.0 / RATE_LIMIT_RPS

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


# ---------------------------
# Rate Limiting Decorator
# ---------------------------

_last_request_time = 0.0

def rate_limited_request(func):
    """Decorator to enforce rate limit and retry on 429 responses."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        global _last_request_time
        elapsed = time.time() - _last_request_time
        if elapsed < MIN_REQUEST_INTERVAL:
            time.sleep(MIN_REQUEST_INTERVAL - elapsed)

        for attempt in range(5):  # up to 5 retries
            response = func(*args, **kwargs)
            if response.status_code != 429:
                _last_request_time = time.time()
                return response

            retry_after = int(response.headers.get("Retry-After", 2))
            logging.warning(f"Rate limited (HTTP 429). Retrying after {retry_after}s...")
            time.sleep(retry_after)

        raise RuntimeError("Exceeded retry limit after repeated 429 responses.")
    return wrapper


@rate_limited_request
def safe_get(url, **kwargs):
    """Wrapper around requests.get with rate limiting and retries."""
    return requests.get(url, **kwargs)


# ---------------------------
# Utilities
# ---------------------------

def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def download_and_extract_zip(url: str, accession: str, output_dir: str):
    """Download a genome ZIP and extract FASTA files."""
    r = safe_get(url, stream=True, headers={"accept": "application/zip"})
    r.raise_for_status()

    with tempfile.NamedTemporaryFile(delete=False, suffix=".zip") as tmp:
        for chunk in r.iter_content(chunk_size=8192):
            tmp.write(chunk)
        tmp_path = tmp.name

    with zipfile.ZipFile(tmp_path, "r") as zf:
        for member in zf.namelist():
            if member.endswith((".fna", ".faa", ".fasta")):
                zf.extract(member, output_dir)
                if member.endswith(".fna"):
                    shutil.move(
                        os.path.join(output_dir, member),
                        os.path.join(output_dir, f"{accession}.fasta")
                    )
                elif member.endswith(".faa"):
                    shutil.move(
                        os.path.join(output_dir, member),
                        os.path.join(output_dir, f"{accession}_proteins.fasta")
                    )

    os.remove(tmp_path)


def translate_6_frames(fasta_path: str, output_path: str):
    """
    Translate nucleotide sequences in fasta_path into 6-frame protein sequences using seqkit.
    Saves the translated sequences to output_path.
    """
    # check that seqkit is installed
    import subprocess
    try:
        subprocess.run(["seqkit", "--version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except (subprocess.CalledProcessError, FileNotFoundError):
        raise RuntimeError("seqkit is not installed or not found in PATH. Please install seqkit to use this feature.")
    # translate in all six frames, append frame info, convert init codon to M and trim trailing X/*
    cmd = ["seqkit", "translate", "-f", "6", "-F", "-M", "--trim", fasta_path]

    out_dir = os.path.dirname(output_path)
    if out_dir:
        ensure_dir(out_dir)

    try:
        with open(output_path, "w") as out_f:
            subprocess.run(cmd, check=True, stdout=out_f, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"seqkit translate failed: {e.stderr.strip()}")


# ---------------------------
# Fetching genomes
# ---------------------------

def fetch_complete_genomes_for_taxon(taxon: str) -> list:
    """Return a list of dicts with assembly info for a taxon at complete genome level."""
    assemblies = []
    next_page_token = None
    page_size = 500  # smaller requests to avoid 504 timeouts

    logging.info(f"Querying NCBI Datasets API for '{taxon}' complete genomes...")

    while True:
        url = f"{API_BASE}/taxon/{requests.utils.quote(taxon)}/dataset_report"
        params = {
            "filters.assembly_level": "complete_genome",
            "page_size": page_size,
        }
        if next_page_token:
            params["page_token"] = next_page_token

        # retry loop for 5xx errors
        for attempt in range(5):
            try:
                r = safe_get(url, params=params, headers={"accept": "application/json"})
                if r.status_code >= 500:
                    raise requests.exceptions.HTTPError(f"{r.status_code} Server Error")
                r.raise_for_status()
                break
            except requests.exceptions.RequestException as e:
                wait_time = 2 ** attempt
                logging.warning(f"Error fetching page ({e}). Retrying in {wait_time}s...")
                time.sleep(wait_time)
        else:
            raise RuntimeError(f"Failed after multiple retries for taxon: {taxon}")

        data = r.json()
        total = data.get("total_count", 0)
        logging.info(f"Found {total} complete genomes for {taxon}")

        for record in data.get("reports", []):
            acc = record.get("accession")
            if acc:
                assemblies.append({"accession": acc})

        next_page_token = data.get("next_page_token")
        if not next_page_token:
            break

        # pause slightly between page requests (stay within rate limit)
        time.sleep(MIN_REQUEST_INTERVAL)

    logging.info(f"Retrieved {len(assemblies)} total assemblies for {taxon}")
    return assemblies


def fetch_random_genomes(taxon: str, n: int, output_dir: str):
    logging.info(f"Fetching list of complete genomes for {taxon}...")
    assemblies = fetch_complete_genomes_for_taxon(taxon)

    if not assemblies:
        raise RuntimeError(f"No complete genomes found for {taxon}")

    chosen = random.sample(assemblies, min(n, len(assemblies)))
    logging.info(f"Selected {len(chosen)} genomes for download.")
    ensure_dir(output_dir)

    for asm in tqdm(chosen, desc="Downloading"):
        acc = asm["accession"]
        url = f"{API_BASE}/accession/{acc}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=PROT_FASTA&hydrated=FULLY_HYDRATED"
        download_and_extract_zip(url, acc, output_dir)

    logging.info(f"Saved genomes in {output_dir}")


def fetch_from_csv(csv_path: str, output_dir: str):
    ensure_dir(output_dir)

    ids = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            embl_id = row.get("RefSeq_ID")
            if embl_id:
                ids.append(embl_id)

    if not ids:
        raise RuntimeError("No valid RefSeq IDs found in CSV.")

    logging.info(f"Fetching {len(ids)} genomes from CSV...")
    for acc in tqdm(ids, desc="Downloading"):
        url = f"{API_BASE}/accession/{acc}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=PROT_FASTA&hydrated=FULLY_HYDRATED"
        try:
            download_and_extract_zip(url, acc, output_dir)
        except Exception as e:
            logging.warning(f"Failed to fetch {acc}: {e}")

    logging.info(f"Saved genomes in {output_dir}")


# ---------------------------
# Main
# ---------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Fetch NCBI genomes (random or CSV-based) and translate them."
    )
    parser.add_argument("pathogen_dir", help="Output directory under data/")
    parser.add_argument("--csv_file", nargs="?", help="CSV file with strain names and IDs")
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    parser.add_argument("--random", dest="random_pathogen", help="Fetch random genomes for pathogen")
    parser.add_argument("--random-num", type=int, default=DEFAULT_RANDOM_NUM)
    parser.add_argument("--output-dir", help="Output directory (default: pathogen_dir/strain_genomes)", default="strain_genomes")

    args = parser.parse_args()
    setup_logging()

    base_dir = os.path.join("data", args.pathogen_dir)
    output_dir = os.path.join(base_dir, args.output_dir)
    ensure_dir(output_dir)

    if args.random_pathogen:
        fetch_random_genomes(args.random_pathogen, args.random_num, output_dir)
    else: 
        if not args.csv_file:
            parser.error("CSV file required when not using --random")
        csv_path = os.path.join(base_dir, args.csv_file)
        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"CSV file not found: {csv_path}")
        fetch_from_csv(csv_path, output_dir)

    shutil.rmtree(os.path.join(output_dir, "ncbi_dataset")) if os.path.exists(os.path.join(output_dir, "ncbi_dataset")) else None  # remove zip extraction dir if exists

    # translate genomes to proteins in 6 frames
    logging.info("Translating genomes to proteins in 6 frames...")
    for file in tqdm(os.listdir(output_dir), desc="Translating"):
        if file.endswith(".fasta") and not file.endswith("_proteins.fasta"):
            fasta_path = os.path.join(output_dir, file)
            protein_output_path = os.path.join(output_dir, f"{os.path.splitext(file)[0]}_6frame_proteins.fasta")
            translate_6_frames(fasta_path, protein_output_path)

    logging.info("Done.")


if __name__ == "__main__":
    main()
