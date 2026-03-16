"""
run_ifnepitope2.py
Runner for IfNePitope2 immunogenicity prediction tool.

Overview:
    - Ensures the required Conda environment and dependencies are present.
    - Applies patches to fix known bugs in the ifnepitope2 package.
    - Runs IfNePitope2 on the provided input FASTA file.
    - Outputs immunogenicity predictions as a CSV file.

Arguments:
    tool_path (Path): Directory containing tools (kept for interface consistency).
    input_fasta (Path): Path to the input FASTA file.
    output_dir (Path): Directory where output will be saved.
    job_type (int, optional): IfNePitope2 job type (default: 1).

Requirements:
    - ext_tools_dependencies.yml (defines Conda environment).
    - pip-installable `ifnepitope2` package inside that environment.
    - Conda available in PATH.

Outputs:
    <output_dir>/<input_fasta_stem>_ifnepitope2.csv   # Immunogenicity prediction results

Author: Nadia
"""
import logging
from pathlib import Path
import subprocess
from tools import common

def run(tool_path: Path, input_fasta: Path, output_dir: Path, job_type: int = 3):
    """
    Run ifnepitope2 prediction tool from within the conda environment.
    Args:
        tool_path (Path): Directory containing tools (not used in this script).
        input_fasta (Path): Path to the input FASTA file.
        output_dir (Path): Directory where output will be saved.
        job_type (int, optional): IfNePitope2 job type (default: 3).
    """

    input_fasta = Path(input_fasta).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{input_fasta.stem}_ifnepitope2.csv"

    cmd = [
        "conda", "run", "-n", common.EXT_TOOLS_ENV_NAME,
        "ifnepitope2",
        "-i", str(input_fasta),
        "-o", str(output_file),
        "-s", "1",              # host human
        "-j", str(job_type),    # job type
        "-d", "2"               # display mode: all peptides
    ]

    logging.info(f"🚀 Running IfNePitope2 on {input_fasta.name}")
    try:
        subprocess.run(cmd, check=True)
        logging.info(f"✅ IfNePitope2 finished: {output_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"❌ IfNePitope2 failed: {e}")
        raise
