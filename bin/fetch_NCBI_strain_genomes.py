#!/usr/bin/env python3
"""
fetch_ncbi_strain_genomes.py

Fetches complete bacterial genomes from NCBI Datasets v2 API.

Dependencies:
    pip install requests tqdm

Features:
    - Fetch random complete genomes for a given taxon.
    - Fetch genomes based on a CSV file containing RefSeq IDs.
    - Automatically downloads and extracts genome FASTA and protein FASTA files.
    - Handles pagination for large datasets.
    - Respects NCBI rate limit (max 5 requests/second).

Usage:
    Fetch random genomes:
        python fetch_ncbi_strain_genomes.py --random "Staphylococcus aureus" \
            --random-num 5 data/staph_aureus

    Fetch genomes from a CSV file:
        python fetch_ncbi_strain_genomes.py data/staph_aureus --csv_file strains.csv
"""

import argparse
import csv
import logging
import os
import random
import shutil
import tempfile
import zipfile
import subprocess
import threading
import time
import requests
from tqdm import tqdm
from requests.adapters import HTTPAdapter, Retry

# ---------------------------
# Configuration
# ---------------------------
API_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome"
DEFAULT_THREADS = 4
DEFAULT_RANDOM_NUM = 5

# Global rate-limit enforcement
RATE_LIMIT_LOCK = threading.Lock()
LAST_REQUEST_TIME = 0
MIN_DELAY = 0.2  # 1 / 5 requests per second


# ---------------------------
# Setup
# ---------------------------
def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def ensure_dir(path: str):
    """Ensure directory exists."""
    os.makedirs(path, exist_ok=True)


# ---------------------------
# Networking Utilities
# ---------------------------
def rate_limited_request(session, method, url, **kwargs):
    """
    Wrapper around session.request() enforcing NCBI's 5 req/s limit.
    Applies a 0.2s delay between consecutive requests.
    """
    global LAST_REQUEST_TIME
    with RATE_LIMIT_LOCK:
        elapsed = time.time() - LAST_REQUEST_TIME
        if elapsed < MIN_DELAY:
            time.sleep(MIN_DELAY - elapsed)
        response = session.request(method, url, **kwargs)
        LAST_REQUEST_TIME = time.time()
    return response


def get_requests_session():
    """Create a persistent session with retries for transient NCBI errors."""
    session = requests.Session()
    retries = Retry(
        total=6,
        backoff_factor=0.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
        raise_on_status=False
    )
    adapter = HTTPAdapter(max_retries=retries)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


# ---------------------------
# Downloading and Extraction
# ---------------------------
def download_and_extract_zip(url: str, accession: str, output_dir: str, session=None):
    """
    Download and extract genome/protein FASTA files.
    Retries automatically on rate limits, server errors, and invalid ZIPs.
    """
    session = session or get_requests_session()

    for attempt in range(5):
        try:
            r = rate_limited_request(session, "GET", url, stream=True, headers={"accept": "application/zip"})
            if r.status_code in (429, 500, 502, 503, 504):
                wait = 2 ** attempt
                logging.warning(f"Server error {r.status_code} for {accession}. Retrying in {wait}s...")
                time.sleep(wait)
                continue
            elif not r.ok:
                r.raise_for_status()

            # Validate content type
            content_type = r.headers.get("Content-Type", "")
            if "zip" not in content_type.lower():
                logging.warning(f"Invalid content for {accession}: {content_type}")
                time.sleep(2 ** attempt)
                continue

            with tempfile.NamedTemporaryFile(delete=False, suffix=".zip") as tmp:
                for chunk in r.iter_content(chunk_size=8192):
                    tmp.write(chunk)
                tmp_path = tmp.name

            # Check if ZIP is valid
            try:
                with zipfile.ZipFile(tmp_path, "r") as zf:
                    for member in zf.namelist():
                        if member.endswith((".fna", ".faa", ".fasta")):
                            zf.extract(member, output_dir)
                            src = os.path.join(output_dir, member)
                            if member.endswith(".fna"):
                                dst = os.path.join(output_dir, f"{accession}.fasta")
                            elif member.endswith(".faa"):
                                dst = os.path.join(output_dir, f"{accession}_proteins.fasta")
                            else:
                                dst = os.path.join(output_dir, member)
                            shutil.move(src, dst)
                os.remove(tmp_path)
                return  # ✅ success
            except zipfile.BadZipFile:
                logging.warning(f"Bad ZIP for {accession} (attempt {attempt+1}). Retrying...")
                os.remove(tmp_path)
                time.sleep(2 ** attempt)

        except Exception as e:
            logging.warning(f"Error fetching {accession}: {e} (attempt {attempt+1})")
            time.sleep(2 ** attempt)

    logging.error(f"Failed to download valid ZIP for {accession} after multiple attempts — skipping.")


# ---------------------------
# Translation
# ---------------------------
def translate_6_frames(fasta_path: str, output_path: str):
    """
    Translate nucleotide sequences into 6-frame protein sequences using seqkit.
    Saves translated sequences to output_path.
    """
    try:
        subprocess.run(["which", "seqkit"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except (subprocess.CalledProcessError, FileNotFoundError):
        raise RuntimeError("seqkit is not installed or not found in PATH.")

    cmd = ["seqkit", "translate", "-f", "6", "-F", "-M", fasta_path]
    ensure_dir(os.path.dirname(output_path))

    try:
        with open(output_path, "w") as out_f:
            subprocess.run(cmd, check=True, stdout=out_f, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"seqkit translate failed: {e.stderr.strip()}")


# ---------------------------
# Fetching Genomes
# ---------------------------
def fetch_complete_genomes_for_taxon(taxon: str) -> list:
    """
    Retrieve all complete genomes for a given taxon using the NCBI Datasets v2 API.
    Handles pagination with adaptive page size reduction and retry/backoff for transient errors.
    """
    assemblies = []
    next_page_token = None
    session = get_requests_session()
    page_size = 1000  # Start large, reduce if timeouts occur

    while True:
        url = f"{API_BASE}/taxon/{requests.utils.quote(taxon)}/dataset_report"
        params = {
            "filters.assembly_level": "complete_genome",
            "filters.assembly_source": "refseq",
            "page_size": page_size,
        }
        if next_page_token:
            params["page_token"] = next_page_token

        for attempt in range(5):
            try:
                r = rate_limited_request(session, "GET", url, params=params, headers={"accept": "application/json"})
            except requests.exceptions.RequestException as e:
                logging.warning(f"Connection error fetching page for {taxon}: {e}. Retrying...")
                time.sleep(5 * (2 ** attempt))
                continue

            # Handle server or rate-limit errors
            if r.status_code in (429, 500, 502, 503, 504):
                wait = 5 * (2 ** attempt)
                logging.warning(f"API error {r.status_code} on {taxon}, retrying in {wait}s...")
                time.sleep(wait)

                # If timeout persists, reduce page size
                if r.status_code == 504 and page_size > 100:
                    old_size = page_size
                    page_size = max(100, page_size // 2)
                    logging.warning(f"Reducing page size from {old_size} to {page_size} due to repeated timeouts.")
                continue

            elif not r.ok:
                r.raise_for_status()

            break  # successful response

        # If still failing after retries, abort
        if not r.ok:
            raise RuntimeError(f"Failed to fetch genome list for {taxon}: {r.status_code}")

        data = r.json()

        # Log total only on first page
        if not next_page_token:
            total = data.get("total_count", 0)
            logging.info(f"Found {total} complete genomes for {taxon}")

        # Collect assemblies
        for record in data.get("reports", []):
            acc = record.get("accession")
            if acc:
                assemblies.append({"accession": acc})

        # Pagination
        next_page_token = data.get("next_page_token")
        if not next_page_token:
            break

    return assemblies


def fetch_random_genomes(taxon: str, n: int, output_dir: str):
    """Fetch and download N random complete genomes for a given taxon."""
    logging.info(f"Fetching list of complete genomes for {taxon}...")
    assemblies = fetch_complete_genomes_for_taxon(taxon)
    if not assemblies:
        raise RuntimeError(f"No complete genomes found for {taxon}")

    chosen = random.sample(assemblies, min(n, len(assemblies)))
    logging.info(f"Selected {len(chosen)} genomes for download.")
    ensure_dir(output_dir)

    session = get_requests_session()
    for asm in tqdm(chosen, desc="Downloading"):
        acc = asm["accession"]
        url = (
            f"{API_BASE}/accession/{acc}/download?"
            "include_annotation_type=GENOME_FASTA&"
            "include_annotation_type=PROT_FASTA&hydrated=FULLY_HYDRATED"
        )
        try:
            download_and_extract_zip(url, acc, output_dir, session=session)
        except requests.exceptions.HTTPError as e:
            logging.warning(f"Skipping {acc}: {e}")

    logging.info(f"Saved genomes in {output_dir}")


def fetch_from_csv(csv_path: str, output_dir: str):
    """Fetch genomes using a CSV file with RefSeq IDs."""
    ensure_dir(output_dir)

    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        ids = [row.get("RefSeq_ID") for row in reader if row.get("RefSeq_ID")]

    if not ids:
        raise RuntimeError("No valid RefSeq IDs found in CSV.")

    session = get_requests_session()
    logging.info(f"Fetching {len(ids)} genomes from CSV...")
    for acc in tqdm(ids, desc="Downloading"):
        url = (
            f"{API_BASE}/accession/{acc}/download?"
            "include_annotation_type=GENOME_FASTA&"
            "include_annotation_type=PROT_FASTA&hydrated=FULLY_HYDRATED"
        )
        try:
            download_and_extract_zip(url, acc, output_dir, session=session)
        except Exception as e:
            logging.warning(f"Failed to fetch {acc}: {e}")

    logging.info(f"Saved genomes in {output_dir}")


# ---------------------------
# Main
# ---------------------------
def main():
    parser = argparse.ArgumentParser(description="Fetch NCBI genomes (random or CSV-based) and translate them.")
    parser.add_argument("pathogen_dir", help="Output directory under data/")
    parser.add_argument("--csv_file", nargs="?", help="CSV file with strain names and IDs")
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    parser.add_argument("--random", dest="random_pathogen", help="Fetch random genomes for pathogen")
    parser.add_argument("--random-num", type=int, default=DEFAULT_RANDOM_NUM)
    parser.add_argument("--output-dir", help="Output directory (default: pathogen_dir/strain_genomes)",
                        default="strain_genomes")

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

    ncbi_dataset_dir = os.path.join(output_dir, "ncbi_dataset")
    if os.path.exists(ncbi_dataset_dir):
        shutil.rmtree(ncbi_dataset_dir)

    # Translate genomes to proteins
    logging.info("Translating genomes to proteins in 6 frames...")
    for file in tqdm(os.listdir(output_dir), desc="Translating"):
        if file.endswith(".fasta") and not file.endswith("_proteins.fasta"):
            fasta_path = os.path.join(output_dir, file)
            protein_output_path = os.path.join(output_dir, f"{os.path.splitext(file)[0]}_6frame_proteins.fasta")
            translate_6_frames(fasta_path, protein_output_path)

    logging.info("Done.")


if __name__ == "__main__":
    main()
