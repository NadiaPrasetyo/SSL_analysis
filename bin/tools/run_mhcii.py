"""
run_mhcii.py
Utility to execute NetMHCIIpan-4.3 epitope prediction tool on a given FASTA input file.

Overview:
    - Runs the external NetMHCIIpan-4.3 tool using subprocess.
    - Uses a predefined set of common MHCII alleles.
    - Organizes output files into a dedicated 'mhcii' subdirectory within the specified output directory.
    - Handles process execution and error reporting.

Arguments:
    fasta_file (str or Path): Path to the input FASTA file containing peptide/protein sequences.
    tool_path (str or Path): Path to the NetMHCIIpan-4.3 executable.
    output_dir (str or Path): Directory where output files will be stored. Results are placed in a 'mhcii' subfolder.

Requirements:
    - NetMHCIIpan-4.3 must be installed and executable at `tool_path`.
    - Python packages: subprocess, pathlib.

Usage Example:
    run("input/protein.fasta",
        "/projects/.../NetMHCIIpan-4.3/Linux_x86_64/bin/netMHCIIpan",
        "results/epitopes")

Outputs:
    - Writes tool output to `output_dir/mhcii/<basename>.out`
    - Writes Excel output to `output_dir/mhcii/<basename>.xls`
    - Prints success or error messages accordingly.

Author: Nadia
"""

import subprocess
from pathlib import Path
import logging


MHCII_DEFAULT = [
    "DRB1_0301", "DRB1_0701", "DRB1_1501",
    "DRB3_0101", "DRB3_0202", "DRB4_0101", "DRB5_0101"
]


def run(fasta_file, tool_path, output_dir):
    """
    Run NetMHCIIpan-4.3 on a given FASTA file.
    Args:
        fasta_file : Path to the input FASTA file containing peptide/protein sequences.
        tool_path : Path to the NetMHCIIpan-4.3 executable.
        output_dir : Directory where output files will be stored. Results are placed in a 'mhcii' subfolder.
    """
    # Default parameters
    rank_threshold = 1.0
    weak_threshold = 5.0
    allele_string = ",".join(MHCII_DEFAULT)

    fasta_file = Path(fasta_file)
    tool_path = Path(tool_path)
    output_dir = Path(output_dir) / "mhcii"
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / f"{fasta_file.stem}_mhcii"

    cmd = [
        str(tool_path),
        "-a", allele_string,
        "-f", str(fasta_file),
        "-rankS", str(rank_threshold),
        "-rankW", str(weak_threshold),
        "-xls",
        "-xlsfile", str(output_file.with_suffix(".xls")),
    ]

    logging.info(f"Running NetMHCIIpan-4.3 on {fasta_file} with {len(MHCII_DEFAULT)} alleles...")
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
        # Save stdout into .out file
        output_file.write_text(result.stdout)

        logging.info(f"✅ NetMHCIIpan-4.3 completed successfully.")
        logging.info(f"Results saved in: {output_file} and {output_file.with_suffix('.xls')}")
    except subprocess.CalledProcessError as e:
        logging.info("❌ NetMHCIIpan-4.3 failed.")
        logging.info("Command:", " ".join(cmd))
        logging.info("Stdout:", e.stdout)
        logging.info("Stderr:", e.stderr)
        raise
