import argparse
from Bio import SeqIO
import csv
import sys

def validate_seq(seq):
    valid_bases = ['A', 'C', 'G', 'T']
    for base in seq:
        if base not in valid_bases:
            return False
    return True

def main(input_fasta, output_file, unique_seq=False):
    # >11849|ERR163737|Staphylococcus_aureus|USA|2016|unknown||8|807188_10863:1-681 + SAUR0420
    # ATGAAATTTAAAGCGATAGCAAAAGCAAGTTTAGCATTGGGAATGTTAGCAACAGGTGTA
    # ATTACATCGAATGTACAATCAGTACAAGCGAAAGCAGAAGTTAAACAACAAAGTGAATCA
    # GAGTTAAAACACTATTATAATAAACCAATTTTAGAGCGTAAAAATGTGACTGGATTTAAA
    # TATACTGATGAGGGTAAACACTATTTAGAAGTCACAGTAGGGCAACAGCATTCTCGAATC
    # ACTTTACTTGGATCTGATAAAGATAAATTTAAAGACGGAGAAAACTCAAATATAGATGTG
    # TTTATCCTTAGAGAAGGTGACAGTAGACAAGCAACAAATTACTCAATTGGTGGCGTTACA
    # AAATCAAATAGTGTGCAGTATATTGATTATATCAATACGCCAATTTTAGAAATCAAGAAA
    # GATAATGAAGATGTACTTAAAGATTTTTACTACATTTCAAAAGAAGACATCTCATTAAAA
    # GAACTTGATTATAGATTAAGAGAACGTGCGATTAAACAACACGGCTTGTATTCAAATGGT
    # CTTAAACAAGGTCAAATTACAATTACAATGAATGATGGCACAACACATACAATCGATTTA
    # AGTCAAAAACTTGAAAAAGAACGTATGGGTGAGTCAATCGACGGCACTAAGATTAATAAA
    # ATTCTAGTAGAAATGAAATAA
    # >ids|isolate|species|country|year|disease|epidemiology|ST (MLST)|seqbin id + position

    csv_data = {}

    with open(input_fasta, "r") as f:
        records = SeqIO.parse(f, "fasta")
        loci = {}
        for record in records:
            isolate = record.description.split("|")[1]
            country = record.description.split("|")[3]
            year = record.description.split("|")[4]
            st = f"ST{record.description.split("|")[7]}"
            location = record.description.split("|")[8]
            locus = location.split("+")[1].strip() if len(location.split("+")) > 1 else "unknown"

            match locus:
                case "SAUR0421":
                    locus = "SSL3"
                case "SAUR0426":
                    locus = "SSL7"
                case "SAUR0420":
                    locus = "SSL11"

            new_header = f">{locus}|{isolate}|{country}|{year}|{st}"
            loci[new_header] = record.seq
            csv_data[isolate] = [locus, country, year, st, record.seq]

    with open(f"{output_file}.csv", "w") as f:
        writer = csv.writer(f)
        writer.writerow(["isolate", "locus", "country", "year", "ST", "sequence"])
        for isolate, data in csv_data.items():
            writer.writerow([isolate] + data)

    with open(f"{output_file}.fasta", "w") as f:
        for header, seq in loci.items():
            f.write(f"{header}\n{seq}\n")

    print(f"Created {output_file}.fasta and {output_file}.csv")


    if unique_seq: # only output unique sequences for each sequence
        unique_seq_file = f"{output_file}_unique.fasta"
        seen_sequences = set()
        with open(unique_seq_file, "w") as f:
            for header, seq in loci.items():
                if seq not in seen_sequences and validate_seq(seq):
                    f.write(f"{header}\n{seq}\n")
                    seen_sequences.add(seq)
        print(f"Created {unique_seq_file}")



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-fasta", required=True, help="Input FASTA file that contains all loci")
    parser.add_argument("-o", "--output", default="pubMLST_summary", help="Output FASTA and CSV file (default: pubMLST_summary)")
    parser.add_argument("--unique-seq", action="store_true", help="output unique sequences to a separate file (default: False)")
    args = parser.parse_args()
    main(args.input_fasta, args.output, args.unique_seq)