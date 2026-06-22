import argparse
import csv
import sys
from Bio import SeqIO

def main(input_file, output_csv):
    # >3532|Saitama9|Staphylococcus_aureus|Japan|2012|1558
    records = SeqIO.parse(input_file, "fasta")

    st_counts = {}
    for record in records:
        st = record.id.split("|")[-1]
        if st in st_counts:
            st_counts[st] += 1
        else:
            st_counts[st] = 1

    with open(output_csv, "w") as f:
        writer = csv.writer(f)
        for st, count in st_counts.items():
            writer.writerow([st, count])

    print(f"Wrote {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", type=str, help="Input file in FASTA format")
    parser.add_argument("-o", "--output_csv", type=str, help="Output CSV file")
    args = parser.parse_args()

    main(args.input_file, args.output_csv)