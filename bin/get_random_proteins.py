from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import random

def get_random_proteins(fasta_file, num_proteins):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    random.shuffle(records)
    return records[:num_proteins]

def main(fasta_file, output_file, num_proteins):
    records = get_random_proteins(fasta_file, num_proteins)
    SeqIO.write(records, output_file, "fasta")

    print(f"Selected {len(records)} proteins from {fasta_file}, saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", type=str, default="random_proteins.fasta", help="Output FASTA file (default: random_proteins.fasta)")
    parser.add_argument("--num", type=int, default=100, help="Number of proteins to select (default: 100)")
    args = parser.parse_args()
    main(args.input, args.output, args.num)