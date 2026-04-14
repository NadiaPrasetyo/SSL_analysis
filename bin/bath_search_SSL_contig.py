import sys
import os
import argparse
import logging
from pathlib import Path
import tempfile
import subprocess

def setup_logging(verbose: bool, output_dir: Path):
    """Configure logging to console and file."""
    log_level = logging.DEBUG if verbose else logging.INFO
    log_file = output_dir / "fetch_pubMLST_contigs.log"
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file, mode='a') if verbose else logging.NullHandler(),
            logging.StreamHandler(sys.stdout)
        ]
    )
    if verbose:
        logging.info("Logging initialized. Log file: %s", log_file)


def check_tool(tool_path):
    if not os.path.exists(tool_path):
        logging.error(f"Tool not found at {tool_path}")
        sys.exit(1)

    if not os.access(tool_path, os.X_OK):
        logging.error(f"Tool at {tool_path} is not executable")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Batch search SSL contig")
    parser.add_argument("-i", "--input", type=str, required=True, help="Input FASTA file")
    parser.add_argument("-t", "--target", type=str, required=True, help="Target FASTA file")
    parser.add_argument("-o", "--output", type=str, default="ssl_contig.txt", help="Output file name (default: ssl_contig.txt)")
    parser.add_argument("--tool-path", type=str, required=True, help="Path to the tool root directory")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    args = parser.parse_args()

    # Set up logging
    output_dir = Path(args.output).parent
    
    setup_logging(args.verbose, output_dir)
    tool_path = args.tool_path
    bathsearch_path = os.path.join(tool_path, "bathsearch")
    bathbuild_path = os.path.join(tool_path, "bathbuild")
    check_tool(bathbuild_path)
    check_tool(bathsearch_path)

    # Run the bath tool to build a hmm profile from the input file
    temp_dir = tempfile.mkdtemp()
    # Use only the stem of the input filename to avoid nested subdirs in temp
    input_stem = Path(args.input).stem
    bhmm_file = os.path.join(temp_dir, f"{input_stem}.bhmm")
    
    #    % bathbuild --unali three_seqs.bhmm three_seqs.fa
    cmd = [
        str(bathbuild_path),
        "--unali",
        bhmm_file,
        args.input
    ]

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logging.error("Error running bathbuild command: %s", e)
        sys.exit(1)

    # run bath search
    #    % bathsearch -o PTH2.out PTH2.bhmm target-PTH2.fa

    cmd = [
        str(bathsearch_path),
        "-o", f"{args.output}.out",
        "--tblout", f"{args.output}.tbl", 
        bhmm_file,
        args.target
    ]

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logging.error("Error running bathsearch command: %s", e)
        sys.exit(1)

    logging.info("Batch search SSL contig completed successfully.")


if __name__ == "__main__":
    main()
