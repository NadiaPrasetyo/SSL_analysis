import subprocess
import tempfile
import argparse
from pathlib import Path
from Bio import SeqIO

def align_mmseqs2(input_fasta, output_dir):
    """
    Runs MMseqs2 easy-search and saves alignment results.
    """

    raw_result = output_dir/ f"{input_fasta}_alignment.tsv"

    with tempfile.TemporaryDirectory() as tmpdir:
        output_fields = [
            "query", "target", "pident", "nident", "alnlen",
            "evalue", "bits", "mismatch", "qcov", "tcov", "tstart", "tend", "taln"
        ]

        db = f"{tmpdir}/db"
        result_db = f"{tmpdir}/result_db"
        result_aln = f"{tmpdir}/result_aln"

        cmd = [
            "mmseqs", "createdb", str(input_fasta), db,
        ]

        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"[✗] MMseqs2 create db failed for {input_fasta}: {e}")
            return None

        cmd = [
            "mmseqs", "align", db, db,
            result_db, result_aln
        ]

        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"[✗] MMseqs2 failed for {input_fasta}: {e}")
            return None
        
        cmd = [
            "mmseqs", "convertalis", result_db, result_db, result_aln, str(raw_result),
            "--format-mode", "4", "--format-output", ",".join(output_fields)
        ]

        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"[✗] MMseqs2 convertalis failed for {input_fasta}: {e}")
            return None

        # ✅ Check for empty or header-only TSV
        has_hits = False
        if raw_result.exists():
            with open(raw_result) as f:
                lines = [ln.strip() for ln in f if ln.strip()]
                if len(lines) > 1:  # more than just a header
                    has_hits = True

        if not has_hits:
            print(f"[✗] No hits found for {input_fasta}")
            return None

        return raw_result
    
    return None

def main():
    parser = argparse.ArgumentParser(description="Align sequences using MMseqs2.")
    parser.add_argument("input_fasta", help="Input FASTA file")
    parser.add_argument("output_dir", help="Output directory")
    args = parser.parse_args()

    align_mmseqs2(args.input_fasta, Path(args.output_dir))

    print(f"[✓] Alignment complete for {args.input_fasta}")


if __name__ == "__main__":
    main()

