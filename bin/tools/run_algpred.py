"""run_algpred.py

Runner for AlgPred2.0 allergenicity prediction tool.

Overview:
    - Ensures the required Conda environment and dependencies are present.
    - Applies patches to fix known bugs in the algpred2 package.
    - Runs AlgPred2.0 on the provided input FASTA file.
    - Outputs allergenicity predictions as a CSV file.

Arguments:
    tool_path (Path): Directory containing tools (kept for interface consistency).
    input_fasta (Path): Path to the input FASTA file.
    output_dir (Path): Directory where output will be saved.
    batch_size (int, optional): Unused, present for interface compatibility.

Requirements:
    - ext_tools_dependencies.yml (defines Conda environment).
    - pip-installable `algpred2` package inside that environment.
    - Conda available in PATH.

Outputs:
    <output_dir>/<input_fasta_stem>_algpred.csv   # Allergenicity prediction results

Author: Nadia
"""

import subprocess
import logging
from pathlib import Path
from tools import common
import py_compile

def patch_algpred_bugs():
    """
    Patch known bugs in algpred2.py:
    - Fix df3.concat misuse
    - Fix str.split error when unpacking 'Name' column
    """
    import sys

    logging.info("🩹 Checking AlgPred2.0 for known bugs...")

    try:
        env_prefix = subprocess.run(
            ["conda", "run", "-n", common.EXT_TOOLS_ENV_NAME, "python", "-c", "import sys; print(sys.prefix)"],
            capture_output=True, check=True, text=True
        ).stdout.strip()

        algpred_path = Path(env_prefix) / "lib/python3.10/site-packages/algpred2/python_scripts/algpred2.py"

        if not algpred_path.exists():
            logging.warning(f"⚠️ Cannot patch: {algpred_path} not found.")
            return

        with open(algpred_path, "r") as f:
            lines = f.readlines()

        new_lines = []
        patched_concat = False
        patched_split = False
        for i, line in enumerate(lines):
            if "df3.concat(" in line:
                leading_ws = line[:len(line) - len(line.lstrip())]
                new_line = f"{leading_ws}df3 = pd.concat([df3, df2.loc[df2.Subject==i][0:5]], axis=0).reset_index(drop=True)\n"
                new_lines.append(new_line)
                patched_concat = True
            elif "df1[['Seq','Hits']] = df1.Name.str.split(" in line:
                leading_ws = line[:len(line) - len(line.lstrip())]
                new_lines.append(f"{leading_ws}split_cols = df1.Name.str.split(\"(\", n=1, expand=True)\n")
                new_lines.append(f"{leading_ws}split_cols.columns = ['Seq', 'Hits']\n")
                new_lines.append(f"{leading_ws}df1 = pd.concat([df1, split_cols], axis=1)\n")
                patched_split = True
            else:
                new_lines.append(line)

        if patched_concat or patched_split:
            logging.info("🔧 Applying patches to algpred2.py...")
            with open(algpred_path, "w") as f:
                f.writelines(new_lines)
            logging.info("✅ Patch applied.")
        else:
            logging.info("✅ No patches needed; bugs not present.")

        # Recompile
        pycache_dir = algpred_path.parent / "__pycache__"
        for pyc_file in pycache_dir.glob("algpred2*.pyc"):
            try:
                pyc_file.unlink()
                logging.info(f"🧹 Removed stale bytecode: {pyc_file.name}")
            except Exception as e:
                logging.warning(f"⚠️ Could not remove bytecode {pyc_file}: {e}")

        try:
            py_compile.compile(str(algpred_path), cfile=None, doraise=True)
            logging.info("🔁 Recompiled algpred2.py successfully.")
        except Exception as e:
            logging.warning(f"⚠️ Bytecode recompilation failed: {e}")

    except Exception as e:
        logging.warning(f"⚠️ Failed to patch AlgPred2.0: {e}")

def run(tool_path: Path, input_fasta: Path, output_dir: Path, batch_size: int = 0):
    """
    Main runner function compatible with pipeline
    Args:
        tool_path (unused): directory containing tools kept for interface consistency
        input_fasta: input FASTA file path
        output_dir: base output directory (tool-specific subdir will be created)
        batch_size (unused): present for interface compatibility
    """
    patch_algpred_bugs()

    input_fasta = Path(input_fasta).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{input_fasta.stem}_algpred.csv"

    cmd = [
        "conda", "run", "-n", common.EXT_TOOLS_ENV_NAME,
        "algpred2",
        "-i", str(input_fasta),
        "-o", str(output_file),
        "-m", "2",  # hybrid-based model
        "-d", "2"   # display mode 2: all peptides (not just allergens)
    ]

    logging.info(f"🚀 Running AlgPred2.0 on {input_fasta.name}")
    try:
        subprocess.run(cmd, check=True)
        logging.info(f"✅ AlgPred2.0 finished: {output_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"❌ AlgPred2.0 failed: {e}")
        raise
