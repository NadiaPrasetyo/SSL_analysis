import requests
import os
import sys
import argparse
from Bio import SeqIO

URL = "https://rest.pubmlst.org/db/"

def fetch_pubMLST_contigs(database, output_dir, date):
    # Fetch the list of isolate records that have genome assemblies
    url = f"{URL}{database}/genomes?return_all=1&added_after={date}"
    response = requests.get(url)
    response.raise_for_status()
    isolate_records = response.json()

    for isolate in isolate_records:
        isolate_id = isolate.split("/")[-1]
        print("Fetching contigs for isolate", isolate_id)
        # Fetch the contigs in FASTA format
        url = f"{URL}{database}/isolates/{isolate_id}/contigs_fasta"
        response = requests.get(url)
        response.raise_for_status()

        # get the contigs into a SeqIO object
        contigs = SeqIO.parse(response.content, "fasta")

        # get other information about the isolate
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch pubMLST contigs")
    parser.add_argument("--database", default="pubmlst_saureus_isolates", help="Database name (e.g.: pubmlst_saureus_isolates)")
    parser.add_argument("--output_dir", default="data/contigs", help="Output directory (default: data/contigs)")
    parser.add_argument("--date", default="2016-01-01", help="Date (YYYY-MM-DD) (default: 2016-01-01)")
    args = parser.parse_args()

    fetch_pubMLST_contigs(args.database, args.output_dir, args.date)
