"""
common.py
Common utility functions and constants for mRNA Target Selection pipeline.

This module provides shared utilities for directory management, file handling, tool validation,
FASTA/JSON conversion, and peptide generation used throughout the mRNA Target Selection workflow.

General Functionality:
    - Ensures output and temporary directories exist and are writable.
    - Locates and validates external tool executables (SignalP, TargetP, TMHMM, IEDB tools, etc.).
    - Provides helper functions for parsing and preparing input/output files.
    - Groups sequences from multiple FASTA files by accession code and merges them.
    - Renames FASTA headers to include only accession codes and strain information.
    - Splits protein FASTA files into peptide FASTA files using a sliding window.
    - Converts FASTA files to text format for downstream processing.
    - Validates peptide sequences based on length and amino acid composition.
    - Prepares output directories for selected tools and cleans up temporary directories.
    - Checks and creates Conda environments for external tools if needed.
    - Retrieves structure files and FASTA files from specified directories.
    - Detects and validates selected antigen tools and IEDB tools.
    - Includes unit tests for key parsing, conversion, and validation functions.

Author: Nadia
"""
from pathlib import Path
import shutil
import tempfile
import logging
import re
from collections import defaultdict
from typing import List
import subprocess
from Bio import SeqIO

def is_valid_peptide(seq: str, min_length: int = 8) -> bool:
    """
    Checks if a given sequence is a valid peptide sequence.
    A valid peptide:
    - Contains only standard amino acids (ACDEFGHIKLMNPQRSTVWY)
    - Is at least `min_length` amino acids long
    """
    return isinstance(seq, str) and len(seq) >= min_length and bool(re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", seq))


def group_cluster_inputs(fasta_files: List[Path], fasta_inputs_dir: Path) -> dict:
    """
    Groups sequences from multiple FASTA files by accession code and merges them.

    Args:
        fasta_files (List[Path]): List of input FASTA files.
        fasta_inputs_dir (Path): Output directory for grouped antigen FASTAs.

    Returns:
        dict: Mapping of accession code to combined FASTA Path.
    """
    grouped_by_accession = defaultdict(list)
    header_pattern = re.compile(
        r"^>[^|]*\|(?P<accession>[A-Z0-9_.-]+)\|[^|]*\|(?P<strain_acc>[A-Z0-9_.-]+)\|.*$"
    )
    # example header formats:
    # >antigen_77|Q5HDD7|Immunoglobulin-binding|HE681097.1|tpos:773380-773815
    # >Q5HDD7|HE681097.1

    for fasta_file in fasta_files:
        with open(fasta_file, "r") as fh:
            lines = fh.readlines()

        i = 0
        while i < len(lines):
            if lines[i].startswith(">"):
                header = lines[i].strip()
                sequence = []
                i += 1
                while i < len(lines) and not lines[i].startswith(">"):
                    sequence.append(lines[i].strip())
                    i += 1

                match = header_pattern.match(header)
                if match:
                    accession = match.group("accession")
                    strain_acc = match.group("strain_acc") or ""  # Default to empty if not present
                    header = f">{accession}|{strain_acc}"  # Reformat header to only have the accession and strain_acc
                    grouped_by_accession[accession].append((header, ''.join(sequence)))
                else:
                    raise ValueError(f"Could not parse accession from header: {header}")
            else:
                i += 1

    output_paths = {}
    for accession, records in grouped_by_accession.items():
        output_path = fasta_inputs_dir / f"{accession}_combined.fasta"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as out_f:
            for header, seq in records:
                out_f.write(f"{header}\n")
                out_f.write(f"{seq}\n")
        output_paths[accession] = output_path

    return output_paths

def rename_fasta_headers(fasta_files: list[Path], tmp_fasta_dir: Path):
    """
    Renames FASTA headers to only include the accession code and strain

    Args:
        fasta_file (Path): Input FASTA file path.
        output_file (Path): Output FASTA file path with renamed headers.
    """
    output_files = []
    for fasta_file in fasta_files:
        output_file = tmp_fasta_dir / fasta_file.name
        with open(fasta_file, "r") as in_f, open(output_file, "w") as out_f:
            for line in in_f:
                # >antigen_77|Q5HDD7|Immunoglobulin-binding|HE681097.1|tpos:773380-773815
                if line.startswith(">"):
                    parts = line.split("|")
                    if len(parts) >= 2:
                        accession = parts[1]
                        strain_acc = parts[3] if len(parts) > 3 else ""  # Default to empty if not present
                        new_header = f">{accession}|{strain_acc}\n"
                        out_f.write(new_header)
                    else:
                        raise ValueError(f"Could not parse accession from header: {line.strip()}")
                else:
                    out_f.write(line)
        output_files.append(output_file)
    
    return output_files

def ensure_writable_dir(path: Path) -> bool:
    """
    Ensure the directory exists and is writable.

    Args:
        path (Path): Directory path to check/create.

    Returns:
        bool: True if the directory exists and is writable, False otherwise.
    """
    try:
        path.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        print(f"❌ Could not create directory {path}: {e}")
        return False

    try:
        with tempfile.TemporaryFile(dir=path):
            pass
    except Exception as e:
        print(f"❌ Directory not writable: {path} ({e})")
        return False

    return True

EXT_TOOLS_ENV_NAME = "external_tools_env"
EXT_TOOLS_ENV_YML = Path("ext_tools_dependencies.yml")

DISCOTOPE_ENV_NAME = "discotope_tools_env"
DISCOTOPE_ENV_YML = Path("discotope_tools_dependencies.yml")

def create_conda_env_if_needed(conda_env_name, conda_env_yml):
    """
    Create Conda environment if it doesn't exist.
    Args:
        conda_env_name (str): Name of the Conda environment.
        conda_env_yml (Path): Path to the Conda environment YAML file.
    Returns:
        None
    """
    logging.info(f"🔍 Checking for Conda environment '{conda_env_name}'...")
    result = subprocess.run(["conda", "env", "list"], capture_output=True, text=True)

    if conda_env_name not in result.stdout:
        logging.info("📦 Environment not found. Creating from YAML...")
        subprocess.run(["conda", "env", "create", "-f", str(conda_env_yml)], check=True)

        # If this is the discotope environment, install PyG packages afterwards
        if conda_env_name == "discotope_tools_env":
            logging.info("📦 Installing PyTorch Geometric packages...")
            subprocess.run([
                "conda", "run", "-n", conda_env_name, "pip", "install",
                "torch-scatter==2.0.9",
                "torch-sparse==0.6.15",
                "torch-geometric==2.1.0.post1",
                "fair-esm==0.5.0"
            ], check=True)

        logging.info("✅ Conda environment created successfully.")
    else:
        logging.info("✅ Conda environment already exists.")

def get_pdb_files(pathogen_path, sequence_dir):
    """
    Get a list of structure files (.pdb, .cif, .cif.gz) from the given sequence directory inside the pathogen path.
    Filters files to include only those with an underscore ("_") in their name.

    Args:
        pathogen_path (Path): Base path to the pathogen.
        sequence_dir (str): Subdirectory under pathogen_path where structure files are located.
    Returns:
        list[Path]: List of filtered structure file paths.
    """
    search_path = pathogen_path / sequence_dir
    extensions = ["*.pdb", "*.cif", "*.cif.gz"]

    structure_files = []
    for ext in extensions:
        structure_files.extend(search_path.glob(ext))

    # Filter files to include only those with an underscore in their name
    filtered_files = [file for file in structure_files if "_" in file.name]

    return filtered_files


def get_fasta_files(base_path: Path, sequence_subdir: str):
    """
    Get all FASTA files in the specified sequence subdirectory.
    Args:
        base_path (Path): Base directory path.
        sequence_subdir (str): Subdirectory name containing FASTA files.
    Returns:
        list: List of Path objects for each FASTA file found.
    """
    seq_dir = base_path / sequence_subdir
    if not seq_dir.exists() or not seq_dir.is_dir():
        print(f"❌ Sequence directory {seq_dir} does not exist or is not a directory.")
        return []
    fasta_files = list(seq_dir.glob("*.fasta"))
    return fasta_files

def prepare_output_dirs(pathogen_path, output_subdir, selected_tools):
    """
    Prepare output directories for the specified pathogen and tools.
    
    Args:
        pathogen_path (Path): Path to the pathogen directory.
        output_subdir (str): Subdirectory name for output files.
        selected_tools (list): List of selected tools to create subdirectories for.
        temp (bool): Whether to create a temporary JSON directory (default: False).
    
    Returns:
        tuple: (output_dir)
    """
    output_dir = pathogen_path / output_subdir

    # Create tool-specific subdirectories under output_dir
    for tool in selected_tools:
        (output_dir / tool.lower()).mkdir(parents=True, exist_ok=True)

    return output_dir

def cleanup_temp(temp_dirs):
    """
    Clean up temporary directories created during processing.
    Args:
        temp_dirs (list): List of temporary Path objects to clean up.
    """
    for temp_dir in temp_dirs:
        if temp_dir.exists():
            try:
                logging.info(f"🗑️ Cleaning up temporary directory: {temp_dir}")
                shutil.rmtree(temp_dir)
            except Exception as e:
                logging.error(f"❌ Failed to delete {temp_dir}: {e}")
        else:
            logging.debug(f"⚠️ Temporary directory does not exist: {temp_dir}")


def check_antigen_tools(tools: list[str], tool_root: Path) -> dict:
    """
    Detects and validates selected antigen tools.
    
    Args:
        tools (list): List of requested tools (e.g., ["SIGNALP", "TARGETP", "CLUSTER", "DEEPLOC", "ALGPRED"])
        tool_root (Path): Root directory where SignalP/TargetP are expected

    Returns:
        dict: Mapping of tool names to either their executable path or placeholder string
    """
    tool_paths = {}

    if "SIGNALP" in tools:
        if not tool_root.exists():
            raise FileNotFoundError(f"Tool root directory {tool_root} does not exist.")
        if not tool_root.is_dir():
            raise NotADirectoryError(f"Tool root {tool_root} is not a directory.")  
        signalp_path = tool_root / "signalp-5.0b" / "bin" / "signalp"
        if not signalp_path.exists():
            raise FileNotFoundError(
                f"SignalP not found at {signalp_path}. "
                "Install via https://services.healthtech.dtu.dk/services/SignalP-5.0/"
            )
        tool_paths["SIGNALP"] = signalp_path

    if "TARGETP" in tools:
        if not tool_root.exists():
            raise FileNotFoundError(f"Tool root directory {tool_root} does not exist.")
        if not tool_root.is_dir():
            raise NotADirectoryError(f"Tool root {tool_root} is not a directory.")
        targetp_path = tool_root / "targetp-2.0" / "bin" / "targetp"
        if not targetp_path.exists():
            raise FileNotFoundError(
                f"TargetP not found at {targetp_path}. "
                "Install via https://services.healthtech.dtu.dk/services/TargetP-2.0/"
            )
        tool_paths["TARGETP"] = targetp_path

    if "CLUSTER" in tools:
        try:
            subprocess.run(["mmseqs", "version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            tool_paths["CLUSTER"] = "mmseqs"  # placeholder
        except Exception:
            raise FileNotFoundError("MMseqs2 not found in PATH. Install via: conda install -c bioconda mmseqs2")

    if "ALGPRED" in tools:
        algpred_path = Path("ext_tools_dependencies.yml")
        if not algpred_path.exists():
            logging.warning(
                f"⚠️ Algpred2 dependencies file not found at {algpred_path}. "
                "Please get the dependency file from GitHub."
            )
        tool_paths["ALGPRED"] = "algpred2"  # placeholder

    if "DEEPLOC" in tools:
        if not tool_root.exists():
            raise FileNotFoundError(f"Tool root directory {tool_root} does not exist.")
        if not tool_root.is_dir():
            raise NotADirectoryError(f"Tool root {tool_root} is not a directory.")  
        deeploc_path = tool_root / "deeplocpro"
        if not deeploc_path.exists():
            raise FileNotFoundError(
                f"Deeplocpro not found at {deeploc_path}. "
                "Install via https://github.com/Jaimomar99/deeplocpro"
            )
        tool_paths["DEEPLOC"] = deeploc_path

    if "IFNEPITOPE2":
        algpred_path = Path("ext_tools_dependencies.yml")
        if not algpred_path.exists():
            logging.warning(
                f"⚠️ Dependencies file not found at {algpred_path}. "
                "Please get the dependency file from GitHub."
            )
        tool_paths["IFNEPITOPE2"] = "ifnepitope2"  # placeholder

    if "DEEPTMHMM":
        # check that pybiolib is installed
        try:
            import biolib
        except ImportError:
            raise ImportError(
                "PyBioLib is not installed. Please install it via: pip install pybiolib"
            )
        tool_paths["DEEPTMHMM"] = "deetmhmm"  # placeholder

    if "MAFFT_RATE4SITE" or "MAFFT" in tools:
        try:
            algpred_path = Path("ext_tools_dependencies.yml")
            if not algpred_path.exists():
                logging.warning(
                    f"⚠️ Dependencies file not found at {algpred_path}. "
                    "Please get the dependency file from GitHub."
                )
            tool_paths["MAFFT_RATE4SITE"] = "mafft"  # placeholder
            tool_paths["MAFFT"] = "mafft"  # placeholder
        except Exception as e:
            raise FileNotFoundError(f"MAFFT_RATE4SITE not found: {e}")

    if "DNDS" in tools:
        try:
            dnds_path = Path("ext_tools_dependencies.yml")
            if not dnds_path.exists():
                logging.warning(
                    f"⚠️ Dependencies file not found at {dnds_path}. "
                    "Please get the dependency file from GitHub."
                )
            tool_paths["DNDS"] = "dnds"  # placeholder
        except Exception as e:
            raise FileNotFoundError(f"DNDS not found: {e}")

    return tool_paths

def check_iedb_tool(base_path):
    """
    Check for the presence of IEDB epitope prediction tools and return their paths.
    Args:
        base_path (str): Base path where IEDB tools are expected to be located.
    Returns:
        dict: Dictionary mapping tool names to their executable paths.
    """
    base = Path(base_path)
    paths = {
        "BCell": base / "BepiPred3_src"/ "bepipred3_CLI.py",
        "MHCI": base / "netMHCpan-4.2" / "netMHCpan",
        "MHCII": base / "netMHCIIpan-4.3"/ "netMHCIIpan",
        "Ellipro": base / "ElliPro.jar",
        "MixMHC2pred": base / "MixMHC2pred-2.0" / "MixMHC2pred_unix",
        "DSSP": base, # just a place holder, it doesn't require a path
        "ProtLearn": base, # just a place holder, it doesn't require a path
        "DiscoTope": base / "discotope3_web" / "src" / "predict_webserver.py"
    }

    tools = {}
    for key, path in paths.items():
        if path.exists():
            tools[key] = str(path)
    return tools

def remove_invalid_aa(fasta_file: Path) -> Path:
    """
    Remove sequences with nonstandard amino acids from a FASTA file.

    Args:
        fasta_file (Path): Path to the FASTA file.

    Returns:
        Path: Path to the cleaned FASTA file (same as input).
    """
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")

    with open(fasta_file) as fin:
        lines = fin.readlines()

    with open(fasta_file, "w") as fout:
        header = None
        seq = ""
        for line in lines:
            if line.startswith(">"):
                if header and seq:
                    if all(aa in valid_aa for aa in seq):
                        fout.write(f"{header}\n{seq}\n")
                    else:
                        logging.warning(f"⚠️ Skipping {header}: contains nonstandard amino acid: {seq}")
                header = line.strip()
                seq = ""
            else:
                seq += line.strip()

        # Write last sequence
        if header and seq and all(aa in valid_aa for aa in seq):
            fout.write(f"{header}\n{seq}\n")

    return fasta_file


def convert_fasta_to_txt(fasta_files, temp_txt_dir: Path):
    """
    Convert FASTA files to text format by copying their contents to .txt files.
    Args:
        fasta_files (list): List of Path objects for FASTA files.
        temp_txt_dir (Path): Directory where the converted text files will be saved.
    Returns:
        list: List of Path objects for the created text files.
    """
    temp_txt_dir.mkdir(parents=True, exist_ok=True)
    cleaned_fasta_files = [remove_invalid_aa(fasta_file) for fasta_file in fasta_files]
    txt_files = []
    for fasta_file in cleaned_fasta_files:
        txt_file = temp_txt_dir / (fasta_file.stem.replace("_cleaned", "") + ".txt")
        txt_file.write_text(fasta_file.read_text())
        txt_files.append(txt_file)
    return txt_files


def split_protein_fasta_to_peptides(input_fasta, output_dir, peptide_length=15):
    """
    Splits one or more protein FASTA files into peptide FASTA files using a sliding window.

    Args:
        input_fasta (str, Path, or list): Path(s) to protein FASTA file(s).
        output_dir (str or Path): Directory to save peptide FASTA files.
        peptide_length (int): Length of peptides (default=15).

    Returns:
        list[Path]: List of generated peptide FASTA file paths.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Allow input_fasta to be a single path or a list of paths
    if isinstance(input_fasta, (str, Path)):
        input_fasta = [Path(input_fasta)]
    else:
        input_fasta = [Path(f) for f in input_fasta]

    output_files = []

    for fasta_file in input_fasta:
        output_fasta = output_dir / (fasta_file.stem + "_peptides.fasta")
        with open(output_fasta, "w") as out_f:
            for record in SeqIO.parse(fasta_file, "fasta"):
                protein_seq = str(record.seq)
                header = record.id

                for i in range(len(protein_seq) - peptide_length + 1):
                    peptide = protein_seq[i:i + peptide_length]
                    peptide_id = f"{header}_seq{i+1}"
                    out_f.write(f">{peptide_id}\n{peptide}\n")

        print(f"✅ Peptide FASTA written to: {output_fasta}")
        output_files.append(output_fasta)

    return output_files


# Add __main__ with unittest
if __name__ == "__main__":
    """ Main function for unit testing of allele retrieval functionality. """
    import unittest

    class RenameFastaHeadersTests(unittest.TestCase):
        def setUp(self):
            self.test_dir = tempfile.TemporaryDirectory()
            self.input_dir = Path(self.test_dir.name) / "input"
            self.output_dir = Path(self.test_dir.name) / "output"
            self.input_dir.mkdir(parents=True)
            self.output_dir.mkdir(parents=True)

            self.fasta1 = self.input_dir / "strain1.fasta"
            self.fasta2 = self.input_dir / "strain2.fasta"

            self.fasta1.write_text(
                ">antigen_153|Q6GHG2|Small|HE681097.1|tpos:380962-381050\n"
                "MAISQERKNEIIKEYRVHETDTGSPEVQIAVLTAEINAVNEHLRTHKKDHHSRRGLLKMVGRRRHLLNYLRSKDIQRYRELIKSLGIRR\n"
            )

            self.fasta2.write_text(
                ">antigen_153|Q6GHG2|Small|HE681098.1|tpos:380962-381050\n"
                "MAISQERKNEIIKEYRVHETDTGSPEVQIAVLTAEINAVNEHLRTTESTHSRRGLLKMVGRRRHLLNYLRSKDIQRYRELIKSLGIRR\n"
                ">antigen_149|Q99QV7|Putative|HE681097.1|tpos:139239-139462\n"
                "MIEFRQVSKTFNKKKQKIHALKDVSFKVNRNDIFGVIGYSGAGKSTLVRLVNHLEAASSGQVLVDGHDITNY\n"
            )

        def tearDown(self):
            self.test_dir.cleanup()

        def test_rename_fasta_headers(self):
            result_files = rename_fasta_headers(
                fasta_files=[self.fasta1, self.fasta2],
                tmp_fasta_dir=self.output_dir
            )
            expected_files = [
                self.output_dir / "strain1.fasta",
                self.output_dir / "strain2.fasta"
            ]
            self.assertEqual(set(result_files), set(expected_files))

            with open(expected_files[0], "r") as f:
                content = f.read()
                self.assertIn(">Q6GHG2|HE681097.1", content)
                self.assertNotIn("antigen_153", content)

            with open(expected_files[1], "r") as f:
                content = f.read()
                self.assertIn(">Q6GHG2|HE681098.1", content)
                self.assertIn(">Q99QV7|HE681097.1", content)
                self.assertNotIn("antigen_153", content)
                self.assertNotIn("antigen_149", content)

    
    class GroupClusterInputsTests(unittest.TestCase):
        def setUp(self):
            self.test_dir = tempfile.TemporaryDirectory()
            self.input_dir = Path(self.test_dir.name) / "input"
            self.output_dir = Path(self.test_dir.name) / "output"
            self.input_dir.mkdir(parents=True)
            self.output_dir.mkdir(parents=True)

            self.fasta1 = self.input_dir / "strain1.fasta"
            self.fasta2 = self.input_dir / "strain2.fasta"

            self.fasta1.write_text(
                ">antigen_153|Q6GHG2|Small|HE681097.1|tpos:380962-381050\n"
                "MAISQERKNEIIKEYRVHETDTGSPEVQIAVLTAEINAVNEHLRTHKKDHHSRRGLLKMVGRRRHLLNYLRSKDIQRYRELIKSLGIRR\n"
            )

            self.fasta2.write_text(
                ">antigen_153|Q6GHG2|Small|HE681098.1|tpos:380962-381050\n"
                "MAISQERKNEIIKEYRVHETDTGSPEVQIAVLTAEINAVNEHLRTTESTHSRRGLLKMVGRRRHLLNYLRSKDIQRYRELIKSLGIRR\n"
                ">antigen_149|Q99QV7|Putative|HE681097.1|tpos:139239-139462\n"
                "MIEFRQVSKTFNKKKQKIHALKDVSFKVNRNDIFGVIGYSGAGKSTLVRLVNHLEAASSGQVLVDGHDITNY\n"
            )

        def tearDown(self):
            self.test_dir.cleanup()

        def test_group_by_accession_combined_output(self):
            result = group_cluster_inputs(
                fasta_files=[self.fasta1, self.fasta2],
                fasta_inputs_dir=self.output_dir
            )

            expected_files = {
                "Q6GHG2": self.output_dir / "Q6GHG2_combined.fasta",
                "Q99QV7": self.output_dir / "Q99QV7_combined.fasta"
            }

            for accession, path in expected_files.items():
                self.assertIn(accession, result)
                self.assertTrue(path.exists())

            # Check Q6GHG2 file for both strain accessions
            with open(expected_files["Q6GHG2"], "r") as f:
                content = f.read()
                self.assertEqual(content.count(">"), 2)
                self.assertIn(">Q6GHG2|HE681097.1", content)
                self.assertIn(">Q6GHG2|HE681098.1", content)

            # Check Q99QV7 file for correct strain accession
            with open(expected_files["Q99QV7"], "r") as f:
                content = f.read()
                self.assertEqual(content.count(">"), 1)
                self.assertIn(">Q99QV7|HE681097.1", content)


    unittest.main()
