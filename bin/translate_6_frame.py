#!/usr/bin/env python3
"""
translate_6frames.py

Translate nucleotide FASTA files into 6-frame protein sequences using seqkit.
Only translates files that do NOT already have corresponding *_proteins.fasta outputs.
"""

import argparse
import os
import subprocess
from tqdm import tqdm


def ensure_dir(path: str):
    """Create directory if it doesn’t exist."""
    if path and not os.path.exists(path):
        os.makedirs(path, exist_ok=True)


def check_seqkit():
    """Ensure seqkit is installed and accessible."""
    try:
        subprocess.run(["seqkit", "version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        raise RuntimeError("❌ seqkit is not installed or not found in PATH. Please install it first.")


def translate_6_frames(fasta_path: str, output_path: str):
    """Translate nucleotide sequences into 6-frame protein sequences using seqkit."""
    ensure_dir(os.path.dirname(output_path))
    cmd = ["seqkit", "translate", "-f", "6", "-F", "-M", fasta_path]

    try:
        with open(output_path, "w") as out_f:
            subprocess.run(cmd, check=True, stdout=out_f, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"❌ seqkit translate failed for {fasta_path}:\n{e.stderr.strip()}")


def translate_all_in_dir(base_dir: str):
    """Translate all FASTA files in a directory that lack corresponding *_proteins.fasta files."""
    fasta_files = [
        f for f in os.listdir(base_dir)
        if f.endswith(".fasta") and not f.endswith("_proteins.fasta")
    ]

    if not fasta_files:
        print("⚠️  No nucleotide FASTA files found for translation.")
        return

    print(f"🔬 Checking {len(fasta_files)} FASTA files for missing translations...")
    to_translate = []

    for file in fasta_files:
        fasta_path = os.path.join(base_dir, file)
        protein_output_path = os.path.join(base_dir, f"{os.path.splitext(file)[0]}_6frame_proteins.fasta")
        if not os.path.exists(protein_output_path):
            to_translate.append((fasta_path, protein_output_path))

    if not to_translate:
        print("✅ All FASTA files already have corresponding _proteins.fasta outputs.")
        return

    print(f"🚀 Translating {len(to_translate)} file(s) into 6-frame proteins...")
    check_seqkit()
    for fasta_path, protein_output_path in tqdm(to_translate, desc="Translating", unit="file"):
        translate_6_frames(fasta_path, protein_output_path)

    print("✅ Translation complete.")


def main():
    parser = argparse.ArgumentParser(description="Translate all FASTA files without _proteins.fasta counterparts.")
    parser.add_argument("sequence_dir", help="Directory containing FASTA files (absolute or relative path)")
    args = parser.parse_args()

    base_dir = args.sequence_dir
    if not os.path.isdir(base_dir):
        raise FileNotFoundError(f"❌ Directory not found: {base_dir}")

    translate_all_in_dir(base_dir)


if __name__ == "__main__":
    main()

