import subprocess
import tempfile
import argparse
from pathlib import Path
from Bio import SeqIO

def align_mmseqs2(input_fasta: Path, output_dir: Path, mode="protein", target_fasta=None):
    """
    Runs MMseqs2 self-alignment for sequences in a FASTA file
    and writes a TSV alignment file.

    Returns path to TSV if hits found, otherwise None.
    """

    input_fasta = Path(input_fasta)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    raw_result = output_dir / f"{input_fasta.stem}_alignment.tsv"

    output_fields = [
        "query", "target", "pident", "alnlen",
        "evalue", "bits", "qcov", "tcov",
        "tstart", "tend", "taln"
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        query_db = tmpdir / "query_db"
        target_db = tmpdir / "target_db"
        result_db = tmpdir / "result_db"

        try:
            # Create DB
            subprocess.run(
                ["mmseqs", "createdb", str(input_fasta), str(query_db)],
                check=True
            )
            if target_fasta is not None:
                subprocess.run(
                    ["mmseqs", "createdb", str(target_fasta), str(target_db)],
                    check=True
                )

            subprocess.run([
                "mmseqs", "search",
                str(query_db), 
                str(query_db) if target_fasta is None else str(target_db),
                str(result_db),
                str(tmpdir),
                "--search-type", "1" if mode == "protein" else "3",
                "-a"
            ], check=True)

            # Convert to TSV
            subprocess.run(
                [
                    "mmseqs", "convertalis",
                    str(query_db), str(query_db),
                    str(result_db),
                    str(raw_result),
                    "--search-type", "1" if mode == "protein" else "3",
                    "--format-mode", "4",
                    "--format-output", ",".join(output_fields),
                ],
                check=True
            )

        except subprocess.CalledProcessError as e:
            print(f"[✗] MMseqs2 failed for {input_fasta}: {e}")
            return None

    # Check if file contains hits
    if raw_result.exists() and raw_result.stat().st_size > 0:
        return raw_result

    print(f"[✗] No hits found for {input_fasta}")
    return None

def main():
    parser = argparse.ArgumentParser(description="Align sequences using MMseqs2.")
    parser.add_argument("input_fasta", help="Input FASTA file")
    parser.add_argument("output_dir", help="Output directory")
    parser.add_argument("--mode", default="protein", help="Alignment mode (protein or nucleotide)")
    parser.add_argument("--target_fasta", default=None, help="Target database for self-alignment")
    args = parser.parse_args()

    align_mmseqs2(args.input_fasta, Path(args.output_dir), mode=args.mode, target_fasta=args.target_fasta)

    print(f"[✓] Alignment complete for {args.input_fasta}")


if __name__ == "__main__":
    main()

