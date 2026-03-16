"""
run_bcell.py
Runner for B-cell epitope prediction tools.

Overview:
    - Applies the BepiPred-3.0 B-cell epitope prediction algorithm to a given FASTA file.
    - Uses a Conda environment ('external_tools_env') to ensure all dependencies are met.
    - Saves prediction results in a structured directory format for downstream analysis.

Arguments:
    fasta_file (Path): Path to the input FASTA file containing protein sequence(s).
    tool_path (str): Path to the BepiPred-3.0 CLI script (e.g., bepipred3_CLI.py).
    output_dir (Path): Directory where prediction results will be saved.

Requirements:
    - Python packages: subprocess, pathlib, logging.
    - A Conda environment named 'external_tools_env' must be created and configured with the required dependencies.
    - The BepiPred-3.0 tool must be installed and accessible.

Outputs:
    <output_dir>/bcell/<fasta_file_stem>/               # Directory containing prediction results.

Notes:
    - This script uses the default prediction mode ('vt_pred') and thresholds provided by BepiPred-3.0.
    - Logs are generated to provide detailed information about the execution process.
    - Ensure that the Conda environment is activated and accessible before running the script.

Author: Nadia
"""
import subprocess
from pathlib import Path
from tools import common
import logging

def run(fasta_file: Path, tool_path: str, output_dir: Path):
    """
    Run BepiPred-3.0 predictor on a given FASTA file using the specified conda environment.

    Args:
        fasta_file : Path to input FASTA file containing protein sequence(s).
        tool_path : Path to the BepiPred-3.0 installation directory (must contain bepipred3_CLI.py).
        output_dir : Path to the output directory where prediction results will be stored.

    Notes:
    - This script assumes that a Conda environment named 'external_tools_env'
      is already created and contains all dependencies listed in requirements.txt.
    - It uses default prediction mode (vt_pred) and default thresholds.
    """

    # Validate input
    fasta_file = Path(fasta_file)
    output_dir = Path(output_dir)/"bcell"/fasta_file.stem
    tool_path = Path(tool_path)
    # tool path would be like: /path/to/BepiPred3_src/bepipred3_CLI.py

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "conda", "run", "-n", common.EXT_TOOLS_ENV_NAME,
        "python", str(tool_path),
        "-i", str(fasta_file),
        "-o", str(output_dir),
        "-pred", "vt_pred"
    ]

    logging.info(f"[INFO] Running BepiPred-3.0:\n{cmd}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        logging.info("STDOUT:\n%s", result.stdout)
        logging.info("STDERR:\n%s", result.stderr)
        
        logging.info(f"[INFO] BepiPred-3.0 run completed successfully. results saved to {output_dir}")
        logging.info(result.stdout)

    except subprocess.CalledProcessError as e:
        logging.error("BepiPred-3.0 failed with exit code %s", e.returncode)
        logging.error("STDOUT:\n%s", e.stdout)
        logging.error("STDERR:\n%s", e.stderr)
        raise      # re-raise or handle as needed


