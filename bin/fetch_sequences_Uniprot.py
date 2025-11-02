"""
fetch_sequences_Uniprot.py
Command-line tool to fetch protein sequence and metadata from UniProt based on antigen data.

Overview:
    - Loads antigen records (antigen names, gene names, UniProt IDs) from a compiled CSV file.
    - Queries the UniProt API to retrieve full protein information for each antigen.
    - Fetches RefSeq nucleotide or protein sequences for matching entries, if available.
    - Parses and standardizes protein metadata including sequence, organism, domains, and features.
    - Compiles and saves the protein data into a new CSV file for downstream analysis.

Arguments:
    pathogen_directory (str): Subdirectory under `data/` containing pathogen data and antigen CSV.
    pathogen_name (str): Full name of the organism (used in UniProt queries).
    --output (str): Optional output CSV file path for compiled protein data. Defaults to
                    `data/<pathogen>/<organism>_compiled_proteins.csv`.
    --input (str): Optional input CSV file path with antigen data. Defaults to
                   `data/<pathogen>/<organism>_compiled_antigens.csv`.
    --fasta (bool): If specified, outputs protein data in FASTA format instead of CSV.

Requirements:
    - Input CSV file with antigen data present in the specified pathogen directory.
    - Python packages: argparse, csv, requests, os, re, unicodedata, time.

Usage Example:
    python fetch_sequences_Uniprot.py sars_cov_2 "SARS-CoV-2" --output proteins.csv --input antigens.csv

Outputs:
    - A CSV file (or FASTA file if --fasta is specified) containing compiled protein metadata for the specified organism, including:
        - UniProt accession
        - Protein name
        - Protein sequence
        - Organism name
        - Pfam domains
        - RefSeq nucleotide or protein sequences (if available)

Author: Nadia
"""
import csv
import requests
import os
import re
import unicodedata
import argparse
import time

def clean_antigen_name(name):
    """
    Clean and standardize the antigen name.
    Args:
        name (str): Original antigen name.
    Returns:
        str: Cleaned antigen name.
    """
    if not isinstance(name, str):
        return ""
    name = name.strip()
    if name.startswith("[") and name.endswith("]"):
        name = name[1:-1].strip()
    if name.startswith("'") and name.endswith("'"):
        name = name[1:-1].strip()

    name = unicodedata.normalize('NFKD', name).encode('ascii', 'ignore').decode('ascii')
    name = name.lower()
    name = re.sub(r'\(.*?\)', '', name)
    name = re.sub(r'\[.*?\]', '', name)
    name = re.sub(
        r'\b(alpha|beta|gamma|delta|epsilon|zeta|theta|kappa|lambda|mu|nu|xi|omicron|pi|rho|sigma|tau|upsilon|phi|chi|psi|omega)[ -]?',
        '', name, flags=re.IGNORECASE
    )
    name = re.split(r'[,/;]', name)[0]
    name = re.sub(r'\s+', ' ', name).strip()
    name = re.sub(r'^\W+|\W+$', '', name)
    return name


def fetch_uniprot_data(query, retries=3, delay=5):
    """Fetch data from UniProt API with retries.
    Args:
        query (str): Query string for UniProt API.
        retries (int): Number of retry attempts on failure.
        delay (int): Delay in seconds between retries.
    Returns:
        dict: JSON response from UniProt API or None if all retries fail.
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": query,
        "fields": "accession,protein_name,sequence,organism_name,xref_pfam,xref_refseq",
        "format": "json",
        "size": 1  # only fetch the top hit
    }

    for attempt in range(retries):
        try:
            response = requests.get(base_url, params=params)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            print(f"Attempt {attempt + 1} failed: {e}")
            if attempt < retries - 1:
                print(f"Retrying in {delay} seconds...")
                time.sleep(delay)
            else:
                print("All retries failed.")
                return None


def fetch_refseq_nucleotide(refseq_id):
    """Fetch RefSeq sequence and return bare nucleotide/protein sequence without headers."""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    db = "protein" if refseq_id.startswith(("WP_", "YP_", "NP_")) else "nuccore"

    params = {
        "db": db,
        "id": refseq_id,
        "rettype": "fasta",
        "retmode": "text"
    }

    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        # Remove FASTA headers and newlines
        lines = response.text.strip().splitlines()
        seq = "".join(line for line in lines if not line.startswith(">"))
        return seq
    except requests.exceptions.RequestException as e:
        print(f"Failed to fetch RefSeq {refseq_id} from {db}: {e}")
        return ""


def parse_uniprot_response(data):
    """Parse UniProt API response to extract relevant fields.
    Args:
        data (dict): JSON response from UniProt API.
    Returns:
        list: List of dictionaries with parsed protein data.
    """
    results = []
    if not data or "results" not in data:
        return results

    for entry in data.get("results", []):
        accession = entry.get("primaryAccession", "")
        protein_name = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
        sequence = entry.get("sequence", {}).get("value", "")
        organism = entry.get("organism", {}).get("scientificName", "")
        pfam = ";".join([xref["id"] for xref in entry.get("uniProtKBCrossReferences", []) if xref.get("database") == "Pfam"])

        # Only fetch the first RefSeq nucleotide sequence
        refseq_ids = [xref["id"] for xref in entry.get("uniProtKBCrossReferences", []) if xref.get("database") == "RefSeq"]
        nucleotide_sequence = ""
        if refseq_ids:
            nucleotide_sequence = fetch_refseq_nucleotide(refseq_ids[0])

        results.append({
            "uniprot_accession": accession,
            "protein_name": protein_name,
            "sequence": sequence,
            "organism_name": organism,
            "pfam": pfam,
            "nucleotide_sequence": nucleotide_sequence
        })
    return results


def load_antigen_records(file_path):
    """Load antigen records from a CSV file.
    Args:
        file_path (str): Path to the CSV file.
    Returns:
        list: List of dictionaries with antigen data.
    """
    records = []
    with open(file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            antigen = {
                "antigen_name": clean_antigen_name(row.get("antigen_name", "")),
                "gene_name": row.get("gene_name", "").strip(),
                "uniprot_id": row.get("Uniprot_ID", "").strip()
            }
            records.append(antigen)
    return records


def main(pathogen, organism="null", output=None, input=None, fasta=False):
    """Main function to fetch and compile protein data from UniProt.
    Args:
        pathogen (str): Subdirectory under `data/` containing pathogen data.
        organism (str): Full name of the organism for UniProt queries.
        output (str): Output CSV file path for compiled protein data.
        input (str): Input CSV file path with antigen data.
        fasta (bool): If True, output FASTA files instead of CSV.
    """
    pathogen_dir = os.path.join("data", pathogen)
    organism_tag = organism.replace(" ", "_").lower() if organism.lower() != "null" else "null"

    antigens_file = os.path.join(pathogen_dir, input or f"{organism_tag}_compiled_antigens.csv")
    output_file = os.path.join(pathogen_dir, output or f"{organism_tag}_compiled_proteins.csv")

    if not os.path.exists(antigens_file):
        print(f"[FATAL] Antigen file not found: {antigens_file}")
        return

    antigen_records = load_antigen_records(antigens_file)
    protein_data = []
    seen_accessions = set()

    for record in antigen_records:
        antigen_name = record.get("antigen_name")
        gene_name = record.get("gene_name")
        uniprot_id = record.get("uniprot_id")

        print(f"[INFO] Processing antigen: {antigen_name}, Gene: {gene_name}, UniProt ID: {uniprot_id}")

        query = None
        if uniprot_id:
            query = f"accession:{uniprot_id}"
        elif gene_name:
            query = f"gene:{gene_name}"
        elif antigen_name:
            query = f"protein_name:{antigen_name}"

        if organism.lower() != "null" and query:
            query += f" AND organism_name:\"{organism}\""

        if query:
            entry = fetch_uniprot_data(query)
            parsed = parse_uniprot_response(entry)
            for p in parsed:
                if p["uniprot_accession"] not in seen_accessions:
                    protein_data.append(p)
                    seen_accessions.add(p["uniprot_accession"])
                    print(f"  [✓] Match: {p['uniprot_accession']}")
                else:
                    print("  [WARN] Duplicate entry skipped")

    if not protein_data:
        print("[WARN] No matching proteins found.")
        return

    fieldnames = ["uniprot_accession", "protein_name", "sequence", "organism_name", "pfam", "nucleotide_sequence"]
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(protein_data)

    print(f"[DONE] Wrote {len(protein_data)} proteins to: {output_file}")


if __name__ == "__main__":
    """Command-line interface for fetching protein sequences from UniProt."""
    parser = argparse.ArgumentParser(description="Fetch protein sequence and metadata from UniProt based on antigen data.")
    parser.add_argument("pathogen_directory", help="Directory name under data/")
    parser.add_argument("pathogen_name", help='Organism name (e.g., "SARS-CoV-2")')
    parser.add_argument("--output", help="Output CSV file path")
    parser.add_argument("--input", help="Input antigen CSV file path")
    args = parser.parse_args()

    main(args.pathogen_directory, args.pathogen_name, output=args.output, input=args.input)
