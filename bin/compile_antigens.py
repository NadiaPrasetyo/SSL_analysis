"""
compile_antigens.py
Command-line tool to compile antigen data from IEDB and literature sources for a given organism.

Overview:
    - Loads IEDB antigen data from a CSV file.
    - Loads literature/patent antigen data from Excel files.
    - Standardizes and merges antigen data from all sources.
    - Outputs a compiled CSV file containing all antigen information.

Arguments:
    pathogen_directory (str): Subdirectory under `data/` containing pathogen data.
    pathogen_name (str): Full organism name (used for file naming and metadata).

Requirements:
    - IEDB antigen CSV file and literature/patent Excel files present in the specified pathogen directory.
    - Python packages: pandas, os, glob, re, argparse.

Usage Example:
    python compile_antigens.py sars_cov_2 "SARS-CoV-2"

Outputs:
    data/<pathogen_directory>/<organism_tag>_compiled_antigens.csv   # Compiled antigen data for the organism

Author: Nadia
"""

import pandas as pd
import os
from glob import glob
import re
import argparse
import sys
print = lambda *args, **kwargs: __builtins__.print(*args, **kwargs, flush=True)


def extract_uniprot_id(iri):
    """
    Extracts the UniProt ID from a given IRI string.
    This function parses a given IRI (Internationalized Resource Identifier) and extracts 
    the UniProt accession if present.
    Args:
        iri (str): The IRI string containing the UniProt ID.
    Returns:
        str: The extracted UniProt ID, or None if not found.
    """
    if pd.isna(iri):
        return None
    match = re.search(r"(?:UNIPROT:)?([A-Z0-9]+)$", iri)
    if match:
        return f"{match.group(1)}"
    return None

def load_iedb_antigens(file_path):
    """
    Loads antigen data from an IEDB CSV file, standardizes the columns, and extracts UniProt IDs.
    This function reads the IEDB CSV file, renames columns to a standard format, extracts
    UniProt IDs from the `parent_source_antigen_iri` column, and adds metadata columns.
    Args:
        file_path (str): Path to the IEDB CSV file.
    Returns:
        pd.DataFrame: Standardized DataFrame with columns:
            - source_organism
            - host_organisms
            - antigen_name
            - gene_name (set to None, as IEDB does not provide this)
            - Uniprot_ID (extracted from parent_source_antigen_iri)
            - source (set to "IEDB")
    """
    if not os.path.exists(file_path):
        print(f"Warning: IEDB file not found: {file_path}. Skipping.")
        return pd.DataFrame()
    df = pd.read_csv(file_path)

    df_out = df.rename(columns={
        "source_organism_names": "source_organism",
        "host_organism_names": "host_organisms",
        "parent_source_antigen_names": "antigen_name"
    })

    # Extract and format Uniprot_ID
    if "parent_source_antigen_iri" in df.columns:
        df_out["Uniprot_ID"] = df["parent_source_antigen_iri"].apply(extract_uniprot_id)
    else:
        df_out["Uniprot_ID"] = None
    df_out["gene_name"] = None  # IEDB does not provide gene names
    df_out["source"] = "IEDB"
    return df_out[["source_organism", "host_organisms", "antigen_name", "gene_name", "Uniprot_ID", "source"]]

def load_literature_antigens(file_path, source_organism):
    """
    Loads antigen data from a literature/patent Excel file, standardizes the columns, and adds metadata.
    This function reads an Excel file containing antigen data, renames columns to a standard format,
    and adds metadata columns for source organism and host organisms.
    Args:
        file_path (str): Path to the literature/patent Excel file.
        source_organism (str): Name of the source organism (e.g. "staphylococcus aureus").
    Returns:
        pd.DataFrame: Standardized DataFrame with columns:
            - source_organism
            - host_organisms (set to "Homo sapiens")
            - antigen_name
            - gene_name
            - Uniprot_ID (set to None, as literature does not provide this)
            - source (set to "literature")
    """
    if not os.path.exists(file_path):
        print(f"Warning: Literature file not found: {file_path}. Skipping.")
        return pd.DataFrame()
    df = pd.read_excel(file_path)
    df_out = pd.DataFrame()
    df.columns = df.columns.str.strip()  # removes leading/trailing spaces
    if "Name" not in df.columns:
        print(f"Warning: 'Name' column not found in {file_path}. Skipping.")
        return pd.DataFrame()
    df_out["antigen_name"] = df["Name"]
    if "Gene" in df.columns :
        if df["Gene"].astype(str).str.contains(r"(?i)\bspy", regex=True).any():
            df_out["gene_name"] = df["Gene"].str.replace(r"(?i)\bspy(\d+)", r"SPy_\1", regex=True)
        else:
            df_out["gene_name"] = df["Gene"]

    else:
        df_out["gene_name"] = None #Gene is optional in the literature

    df_out["source_organism"] = source_organism
    df_out["host_organisms"] = "Homo sapiens"
    df_out["Uniprot_ID"] = df["Uniprot_ID"] if "Uniprot_ID" in df.columns else None
    df_out["source"] = "literature"
    return df_out[["source_organism", "host_organisms", "antigen_name", "gene_name", "Uniprot_ID", "source"]]


def normalize_name(name):
    """Lowercase, strip punctuation except spaces, normalize spaces."""
    if not name or pd.isna(name):
        return ""
    name = str(name).lower()
    name = re.sub(r"[^a-z0-9\s]", "", name)  # remove punctuation
    name = re.sub(r"\s+", " ", name).strip()
    return name

def expand_names(antigen_name, gene_name):
    """
    Extracts all possible names from antigen_name and gene_name.
    Handles parentheses, 'also known as', commas, and 'and'.
    """
    names = set()

    # From antigen_name
    if pd.notna(antigen_name):
        text = str(antigen_name)

        # 1. Full antigen name
        names.add(normalize_name(text))

        # 2. Extract from parentheses
        paren_content = re.findall(r"\((.*?)\)", text)
        for content in paren_content:
            # Remove "also known as"
            if "also known as" in content.lower():
                content = re.sub(r"also known as", "", content, flags=re.IGNORECASE)
            for part in re.split(r",| and ", content):
                names.add(normalize_name(part))

        # 3. Remove parentheses content, split by comma/and
        no_paren_text = re.sub(r"\(.*?\)", "", text)
        for part in re.split(r",| and ", no_paren_text):
            names.add(normalize_name(part))

    # From gene_name
    if pd.notna(gene_name):
        for part in re.split(r",| and |\s", str(gene_name)):
            names.add(normalize_name(part))

    return {n for n in names if n}

def check_for_duplicates(combined_df, gene_only=False):
    """
    Finds partial and exact duplicates by expanded antigen+gene names or just genes.
    Args:
        combined_df (pd.DataFrame): The DataFrame containing antigen data.
        gene_only (bool): If True, only checks for duplicates based on gene names.
    Returns:
        pd.DataFrame: The DataFrame with duplicates resolved.
    """
    if gene_only:
        name_sets = {
            idx: {normalize_name(row.get('gene_name', None))}
            for idx, row in combined_df.iterrows() if pd.notna(row.get('gene_name', None))
        }
    else:
        name_sets = {
            idx: expand_names(row['antigen_name'], row.get('gene_name', None))
            for idx, row in combined_df.iterrows()
        }

    visited = set()
    duplicates_groups = []

    # Find overlaps only for rows with valid name_sets
    for i in combined_df.index:
        if i in visited or i not in name_sets:
            continue
        group = {i}
        for j in combined_df.index:
            if i == j or j in visited or j not in name_sets:
                continue
            if name_sets[i] & name_sets[j]:  # non-empty intersection
                group.add(j)
        if len(group) > 1:
            duplicates_groups.append(group)
            visited |= group

    # Interactive resolution (multi-keep allowed)
    for group in duplicates_groups:
        print("\n⚠ Possible duplicate entries found:\n")
        for idx in sorted(group):
            print(
                f"[{idx}] {combined_df.loc[idx, 'antigen_name']} "
                f"| Gene: {combined_df.loc[idx, 'gene_name']} "
                f"| Source: {combined_df.loc[idx, 'source']}"
            )

        keep = set()
        while not keep.issubset(group) or not keep:
            try:
                keep_input = input(
                    f"Enter the indices to KEEP from {sorted(group)} (comma-separated, or 'all' to keep all, 'none' to remove all): "
                ).strip().lower()
                if keep_input == "all":
                    keep = group
                elif keep_input == "none":
                    keep = set()
                    break
                else:
                    keep = {int(x.strip()) for x in keep_input.split(",") if x.strip()}
            except ValueError:
                print("Invalid input. Enter valid indices separated by commas, or 'all'/'none'.")

        drop_ids = group - keep
        combined_df = combined_df.drop(index=drop_ids)
        print(f"✅ Kept indices {sorted(keep)}, removed {len(drop_ids)} duplicates.")

    return combined_df

def split_and_expand_entries(df):
    """
    Expands rows so that each antigen_name and gene_name combination is separate.
    Splits on commas, 'and', and parentheses. Cleans names for search purposes.
    Removes (UniProt:____) from antigen names before processing.
    """
    expanded_rows = []
    
    for _, row in df.iterrows():
        antigen_names = set()
        gene_names = set()
        
        # Only process rows where the source is not "IEDB"
        if row["source"] != "IEDB":
            # --- Extract antigen names ---
            if pd.notna(row["antigen_name"]):
                text = str(row["antigen_name"])
                
                # Remove any (UniProt:XXXX) from text
                text = re.sub(r"\(.*?UniProt:[A-Z0-9]+\)", "", text, flags=re.IGNORECASE)
                
                # Add the cleaned full name
                antigen_names.add(text.strip())
                
                # Split inside parentheses (ignore UniProt content)
                for content in re.findall(r"\((.*?)\)", text):
                    if re.search(r"uniprot:[A-Z0-9]+", content, flags=re.IGNORECASE):
                        continue
                    content = re.sub(r"also known as", "", content, flags=re.IGNORECASE)
                    for part in re.split(r",| and ", content):
                        part = part.strip()
                        if part:
                            antigen_names.add(part)
                
                # Split without parentheses
                no_paren = re.sub(r"\(.*?\)", "", text)
                for part in re.split(r",| and ", no_paren):
                    part = part.strip()
                    if part:
                        antigen_names.add(part)
            
            # --- Extract gene names ---
            if pd.notna(row["gene_name"]):
                for part in re.split(r",| and |\s", str(row["gene_name"])):
                    part = part.strip()
                    if part:
                        gene_names.add(part)
            else:
                gene_names.add(None)
        
        # --- Expand to all name/gene combinations ---
        if antigen_names or gene_names:
            for name in antigen_names:
                for gene in gene_names:
                    new_row = row.copy()
                    new_row["antigen_name"] = name
                    new_row["gene_name"] = gene
                    expanded_rows.append(new_row)
        else:
            # If no expansion is needed, keep the row as is
            expanded_rows.append(row)
    
    return pd.DataFrame(expanded_rows).reset_index(drop=True)

def main(short_name, long_name, output_file):
    """
    Main entry point to compile various sources of antigens.
    """
    organism_tag = str(long_name).replace(" ", "_").lower()
    base_path = f"data/{short_name}"

    # Load IEDB
    iedb_file = os.path.join(base_path, f"{organism_tag}_IEDB_antigens.csv")
    iedb_df = load_iedb_antigens(iedb_file)

    # Load literature/patents (XLSX files)
    literature_dfs = []
    for xlsx_file in glob(os.path.join(base_path, "*.xlsx")):
        lit_df = load_literature_antigens(xlsx_file, long_name)
        if not lit_df.empty:
            literature_dfs.append(lit_df)
    literature_df = pd.concat(literature_dfs, ignore_index=True) if literature_dfs else pd.DataFrame()

    # Combine all sources
    combined_df = pd.concat([iedb_df, literature_df], ignore_index=True)
    if combined_df.empty:
        print("No antigen data found. Exiting.")
        return
    
    print(f"Loaded {len(iedb_df)} IEDB antigens and {len(literature_df)} literature/patent antigens: total {len(combined_df)} entries.")

    # Expand multiple names/genes into separate rows
    combined_df = split_and_expand_entries(combined_df)

    # Remove exact duplicate rows
    combined_df = combined_df.drop_duplicates(subset=["source_organism", "host_organisms", "antigen_name", "gene_name", "Uniprot_ID", "source"])

    # Cross-check duplicates using expanded name sets (names + genes)
    combined_df = check_for_duplicates(combined_df, False)

    # Another check of just genes
    combined_df = check_for_duplicates(combined_df, True)

    print(f"After removing duplicates, {len(combined_df)} unique antigen entries remain.")

    # Output
    if output_file is None:
        output_file = os.path.join(base_path, f"{organism_tag}_compiled_antigens.csv")
    else:
        output_file = os.path.join(base_path, output_file)
    combined_df.to_csv(output_file, index=False)
    print(f"Compiled antigen data saved to: {output_file}")

if __name__ == "__main__":
    """
    Entry point for command-line execution.
    Parses command-line arguments for pathogen directory and name, then calls main function."""
    parser = argparse.ArgumentParser(
        description="Compile antigen data from IEDB and literature sources for a given organism.",
        usage="python compile_antigens.py <pathogen_directory> <pathogen_name> [--output <output_file>]"
    )
    parser.add_argument("pathogen_directory", help="Directory name under data/")
    parser.add_argument("pathogen_name", help='Prefix used in filenames (e.g., "staphylococcus aureus")')
    parser.add_argument("--output", help="Path to the output CSV file", default=None)
    args = parser.parse_args()

    if not os.path.exists(f"data/{args.pathogen_directory}"):
        os.makedirs(f"data/{args.pathogen_directory}")
        print(f"Created directory: data/{args.pathogen_directory}")

    main(args.pathogen_directory, args.pathogen_name, args.output)
