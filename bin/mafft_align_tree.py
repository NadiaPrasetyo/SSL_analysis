"""
run_mafft.py
Runner for MAFFT multiple sequence alignment tool.
Overview:
    - Ensures the required Conda environment and dependencies are present.
    - Runs MAFFT (L-INS-i algorithm) on the provided input FASTA file using the external_tools_env Conda environment.
    - Outputs the aligned sequences in Clustal format, along with a guide tree.
    - Restores original headers in the alignment and tree files.
Arguments:
    tool_path (Path): Path to the MAFFT executable (ignored, required by interface).
    input_fasta (Path): Path to the input FASTA file.
    output_dir (Path): Directory where output will be saved.
    batch_size (int, optional): Unused, present for interface compatibility.
Requirements:
    - ext_tools_dependencies.yml (defines Conda environment).
    - MAFFT installed in the specified Conda environment.
    - Conda available in PATH.
Outputs:
    <output_dir>/<input_fasta_stem>_aligned.fasta   # Aligned sequences in Clustal format.
    <output_dir>/<input_fasta_stem>.tree            # Guide tree in Newick format (if generated).
Author: Nadia
"""
import subprocess
from pathlib import Path
from typing import Optional
import logging
import shutil
from tools import common
import re

############################ HELPER FUNCTIONS ############################

def count_sequences(input_fasta: Path) -> int:
    """
    Count the number of sequences in a FASTA file.
    Args:
    - input_fasta: Path to the input FASTA file.
    Returns:
    - int: Number of sequences in the FASTA file.
    """
    with open(input_fasta, 'r') as f:
        return sum(1 for line in f if line.startswith('>'))
    

def rename_tree_headers_from_files(tree_file: Path, fasta_file: Path):
    """
    Replace MAFFT-style tree IDs with the exact FASTA headers.
    Args:
    - tree_file: Path to the Newick tree file.
    - fasta_file: Path to the aligned FASTA file.
    """
    # Read tree and FASTA contents
    tree_str = tree_file.read_text()
    fasta_str = fasta_file.read_text()

    # Extract FASTA headers as a list
    fasta_headers = [line[1:].strip() for line in fasta_str.splitlines() if line.startswith(">")]

    # Build a map from MAFFT tree ID (numeric prefix + '_' + header with | and . replaced by _) -> original FASTA header
    # Example: For header "C1CEN2|NZ_CP063829.1"
    # MAFFT ID would be like "2_C1CEN2_NZ_CP063829_1"
    # We'll map "C1CEN2_NZ_CP063829_1" (without numeric prefix) to original header
    header_map = {}
    for header in fasta_headers:
        # Transform FASTA header to MAFFT ID pattern (without numeric prefix)
        mafft_id = header.replace("|", "_").replace(".", "_")
        header_map[mafft_id] = header

    # Regex to find MAFFT IDs in tree:
    # IDs look like: (numeric prefix)_<mafft_id>
    # numeric prefix is digits; mafft_id is letters/digits/underscores
    pattern = re.compile(r"\b(\d+_([A-Za-z0-9_]+))\b")

    def replacer(match):
        full_id = match.group(1)    # e.g. "2_C1CEN2_NZ_CP063829_1"
        mafft_id = match.group(2)   # e.g. "C1CEN2_NZ_CP063829_1"
        return header_map.get(mafft_id, full_id)  # Replace with original FASTA header or keep as is

    updated_tree = pattern.sub(replacer, tree_str)

    # Write updated tree back
    tree_file.write_text(updated_tree)
    
############################ RUN MAFFT FUNCTION ############################

def run(input_fasta: Path, output_dir: Path):
    """
    Runs MAFFT using the external_tools_env conda environment.
    Args:
    - tool_path: Optional[Path] to the MAFFT executable (ignored, required by interface); can be None.
    - input_fasta: Path to the input FASTA file.
    - output_dir: Path to the output directory.
    Raises:
    - RuntimeError: If Conda is not available or if the MAFFT command fails.
    Outputs:
    - <output_dir>/<input_fasta_stem>_aligned.fasta: Aligned sequences in Clustal format.
    - <output_dir>/<input_fasta_stem>.tree: Guide tree in Newick format (if generated).
    """
    # Check if the conda environment exists
    common.create_conda_env_if_needed(common.EXT_TOOLS_ENV_NAME, common.EXT_TOOLS_ENV_YML)

    # Check if the input FASTA file contains more than one sequence
    sequence_count = count_sequences(input_fasta)
    if sequence_count <= 1:
        logging.info(f"ℹ️ Input FASTA file {input_fasta.name} contains only one sequence. Alignment is not required. Exiting gracefully.")
        return

    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / (input_fasta.stem + "_aligned.fasta")

    command = [
        "conda", "run", "-n", common.EXT_TOOLS_ENV_NAME,
        "mafft", 
        "--reorder", # Output order: aligned.
        "--treeout", # Guide tree is output to the input.tree file
        "--amino", # Assume the sequences are amino acid
        str(input_fasta)  # Input FASTA file
    ]

    logging.info(f"🔍 Running MAFFT on {input_fasta}...")

    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)

        # Save the aligned output
        with open(output_file, 'w') as f:
            f.write(result.stdout)

        logging.info(f"✅ MAFFT alignment completed. Output saved to {output_file}")

        # Check if a tree file was generated
        temp_file = input_fasta.with_suffix(input_fasta.suffix + ".tree")  # Input file with .tree suffix
        if not temp_file.exists():
            logging.warning("⚠️ MAFFT did not generate a tree file.")
        
        #move the tree file to the output directory if it exists
        final_output_tree_file = output_dir / f"{input_fasta.stem}.tree"
        shutil.move(temp_file, final_output_tree_file)
        rename_tree_headers_from_files(final_output_tree_file, output_file)

        logging.info("🔁 Restored original headers in alignment and tree.")
        logging.info("✅ MAFFT alignment and tree restoration completed. Output files: "
                    f"{output_file}, {final_output_tree_file}")
        
        cleanup_files = [temp_file]
        for f in cleanup_files:
            if f.exists():
                f.unlink()

    except subprocess.CalledProcessError as e:
        logging.error("❌ Error running MAFFT:")
        logging.error(e.stderr)
        raise

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run MAFFT multiple sequence alignment.")
    parser.add_argument("input_fasta", type=Path, help="Path to the input FASTA file.")
    parser.add_argument("output_dir", type=Path, help="Directory where output will be saved.")

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    run(args.input_fasta, args.output_dir)
