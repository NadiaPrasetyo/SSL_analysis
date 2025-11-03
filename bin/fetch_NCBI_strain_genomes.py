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
import logging
import requests
from tqdm import tqdm

# ---------------------------
# Configuration
# ---------------------------
API_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome"
DEFAULT_THREADS = 4
DEFAULT_RANDOM_NUM = 5

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


# ---------------------------
# Utilities
# ---------------------------

def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


import time
from requests.adapters import HTTPAdapter, Retry

def get_requests_session():
    """Create a persistent session with retries for transient NCBI errors."""
    session = requests.Session()
    retries = Retry(
        total=6,
        backoff_factor=0.5,  # quick exponential backoff
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
        raise_on_status=False
    )
    adapter = HTTPAdapter(max_retries=retries)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session

import time
from requests.adapters import HTTPAdapter, Retry

def get_requests_session():
    """Create a persistent session with retries for transient NCBI errors."""
    session = requests.Session()
    retries = Retry(
        total=6,
        backoff_factor=0.5,  # quick exponential backoff
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
        raise_on_status=False
    )
    adapter = HTTPAdapter(max_retries=retries)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def download_and_extract_zip(url: str, accession: str, output_dir: str, session=None):
    """
    Download and extract genome/protein FASTA files.
    Retries automatically on rate limits, server errors, and invalid ZIPs.
    """
    session = session or get_requests_session()

    for attempt in range(5):
        try:
            r = session.get(url, stream=True, headers={"accept": "application/zip"})
            if r.status_code in (429, 500, 502, 503, 504):
                wait = 2 ** attempt
                logging.warning(f"Server error {r.status_code} for {accession}. Retrying in {wait}s...")
                time.sleep(wait)
                continue
            elif not r.ok:
                r.raise_for_status()

            # Validate Content-Type
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
                continue

        except Exception as e:
            logging.warning(f"Error fetching {accession}: {e} (attempt {attempt+1})")
            time.sleep(2 ** attempt)

    logging.error(f"Failed to download valid ZIP for {accession} after multiple attempts — skipping.")

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
    """
    Return a list of dicts with assembly info for a taxon at complete genome level.
    Handles pagination to retrieve all available assemblies, with retry/backoff for 500/429 errors.
    """
    assemblies = []
    next_page_token = None
    session = get_requests_session()

    while True:
        url = f"{API_BASE}/taxon/{requests.utils.quote(taxon)}/dataset_report"
        params = {
            "filters.assembly_level": "complete_genome",
            "page_size": 1000,  # maximum allowed page size
        }
        if next_page_token:
            params["page_token"] = next_page_token

        # retry loop for transient errors
        for attempt in range(5):
            r = session.get(url, params=params, headers={"accept": "application/json"})
            if r.status_code in (429, 500, 502, 503, 504):
                wait = 2 ** attempt
                logging.warning(f"NCBI API error {r.status_code} on {taxon}, retrying in {wait}s...")
                time.sleep(wait)
                continue
            elif not r.ok:
                r.raise_for_status()
            break

        # If still failing after retries, give up
        if not r.ok:
            raise RuntimeError(f"Failed to fetch genome list for {taxon} after retries: {r.status_code}")

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

    return assemblies


def fetch_random_genomes(taxon: str, n: int, output_dir: str):
    logging.info(f"Fetching list of complete genomes for {taxon}...")
    assemblies = fetch_complete_genomes_for_taxon(taxon)
    if not assemblies:
        raise RuntimeError(f"No complete genomes found for {taxon}")

    chosen = random.sample(assemblies, min(n, len(assemblies)))
    logging.info(f"Selected {len(chosen)} genomes for download.")
    logging.info(f"Chosen assemblies: {chosen}")
    ensure_dir(output_dir)

    session = get_requests_session()
    logging.info(f"Downloading {len(chosen)} random genomes...")
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
    ensure_dir(output_dir)

    ids = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            embl_id = row.get("RefSeq_ID")
            if embl_id:
                ids.append(embl_id)

    if not ids:
        raise RuntimeError("No valid EMBL/GenBank IDs found in CSV.")

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
