import requests
import os
import sys
import argparse
import io
import time
from Bio import SeqIO
import logging
from pathlib import Path
import traceback

def setup_logging(verbose: bool, output_dir: Path):
    """Configure logging to console and file."""
    log_level = logging.DEBUG if verbose else logging.INFO
    log_file = output_dir / "fetch_pubMLST_contigs.log"
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

URL = "https://rest.pubmlst.org/db/"

def fetch_pubMLST_contigs(database, output_dir, date):
    url = f"{URL}{database}/genomes?return_all=1&added_after={date}"
    logging.info(f"Fetching isolate records with URL: {url}")
    response = requests.get(url)
    response.raise_for_status()
    response_json = response.json()
    isolate_records = response_json["isolates"]

    for isolate in isolate_records:
        isolate_id = isolate.split("/")[-1]
        logging.info(f"Fetching contigs for isolate {isolate_id}")

        # Fetch the contigs in FASTA format
        url = f"{URL}{database}/isolates/{isolate_id}/contigs_fasta"
        response = requests.get(url)
        response.raise_for_status()

        # ✅ FIX: Wrap bytes in BytesIO so SeqIO treats it as a stream, not a filename
        fasta_io = io.StringIO(response.content.decode("utf-8"))
        contigs = list(SeqIO.parse(fasta_io, "fasta"))  # ✅ FIX: materialise generator into list

        # Fetch isolate metadata
        logging.info(f"Fetching isolate information for isolate {isolate_id}")
        url = f"{URL}{database}/isolates/{isolate_id}"
        response = requests.get(url)
        response.raise_for_status()
        isolate_info = response.json()

        accession = isolate_info.get("provenance", {}).get("run_accession", "unknown")
        country = isolate_info.get("provenance", {}).get("country", "unknown")
        year = isolate_info.get("provenance", {}).get("year", "unknown")
        if year == "unknown":
            year = isolate_info.get("provenance", {}).get("date_entered", "unknown").split("-")[0]

        # schemes is a list — find the MLST entry by description
        st = "unknown"
        for scheme in isolate_info.get("schemes", []):
            if scheme.get("description") == "MLST":
                st = scheme.get("fields", {}).get("ST", "unknown")
                break

        for record in contigs:
            record.id = f"{isolate_id}|{record.id}|{accession}|{country}|{year}|ST-{st}"
            record.description = ""  # prevent duplicate info in FASTA header

        output_file = os.path.join(output_dir, f"{isolate_id}.fasta")
        with open(output_file, "w") as f:
            SeqIO.write(contigs, f, "fasta")
        logging.info(f"Saved contigs for isolate {isolate_id} to {output_file}")

        time.sleep(1)  # ✅ FIX: time is now imported

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch pubMLST contigs")
    parser.add_argument("--database", default="pubmlst_saureus_isolates")
    parser.add_argument("--output_dir", default="data/contigs")
    parser.add_argument("--date", default="2016-01-01")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    setup_logging(args.verbose, Path(args.output_dir))

    try:
        fetch_pubMLST_contigs(args.database, args.output_dir, args.date)
    except Exception as e:
        logging.error(f"Error: {e}")
        logging.debug(traceback.format_exc())