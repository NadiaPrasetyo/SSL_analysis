import requests
import os
import sys
import argparse
from Bio import SeqIO
import logging

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
    # Fetch the list of isolate records that have genome assemblies
    url = f"{URL}{database}/genomes?return_all=1&added_after={date}"
    logging.info("Fetching isolate records with URL:", url)
    response = requests.get(url)
    response.raise_for_status()
    response_json = response.json()
    isolate_records = response_json["isolates"]

    for isolate in isolate_records:
        isolate_id = isolate.split("/")[-1]
        logging.info("Fetching contigs for isolate", isolate_id)
        # Fetch the contigs in FASTA format
        url = f"{URL}{database}/isolates/{isolate_id}/contigs_fasta"
        response = requests.get(url)
        response.raise_for_status()

        # get the contigs into a SeqIO object
        contigs = SeqIO.parse(response.content, "fasta")

        # get other information about the isolate
        logging.info("Fetching isolate information for isolate", isolate_id)
        url = f"{URL}{database}/isolates/{isolate_id}"
        response = requests.get(url)
        response.raise_for_status()
        isolate_info = response.json()

        # change the headers of the contigs to include the isolate information
        # isolate,locus,country,year,clonal_complex,sequence
        country = isolate_info["country"]
        year = isolate_info["year"]
        clonal_complex = isolate_info["clonal_complex"]


        for record in contigs:
            record.id = f"{isolate_id}|{record.id}|{country}|{year}|{clonal_complex}"

        output_file = os.path.join(output_dir, f"{isolate_id}.fasta")
        with open(output_file, "w") as f:
            SeqIO.write(contigs, f, "fasta")

        logging.info(f"Saved contigs for isolate {isolate_id} to {output_file}")
        # add a time out to prevent overloading the server
        time.sleep(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch pubMLST contigs")
    parser.add_argument("--database", default="pubmlst_saureus_isolates", help="Database name (e.g.: pubmlst_saureus_isolates)")
    parser.add_argument("--output_dir", default="data/contigs", help="Output directory (default: data/contigs)")
    parser.add_argument("--date", default="2016-01-01", help="Date (YYYY-MM-DD) (default: 2016-01-01)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    args = parser.parse_args()

    setup_logging(args.verbose, Path(args.output_dir))
    fetch_pubMLST_contigs(args.database, args.output_dir, args.date)
