import argparse
from Bio import SeqIO
import csv
from collections import defaultdict


def main(input_fasta, nail, output):
    # trim the fasta sequences to between the nail target start and end
    # NAIL:
    #                                                                     target target query query       comp          cell  
    # target                                                    query     start  end    start end   score bias evalue   frac  
    # --------------------------------------------------------- --------- ------ ------ ----- ----- ----- ---- -------- ----- 
    nail_start_end = defaultdict(list)
    with open(nail, "r") as f:
        for line in f:
            row_list = line.split()
            if row_list[0] == "#":
                continue
            else:
                if len(row_list) < 4:
                    print("Not enough columns in NAIL file row:", row_list)
                    continue
                try:
                    nail_start_end[row_list[0]] = [int(row_list[2]), int(row_list[3])]
                except ValueError:
                    if row_list[0] == "target":
                        raise
                    else:
                        print("Invalid literal for int() with base 10:", row_list)
                        raise

    # trim the fasta sequences in the fasta file
    records = list(SeqIO.parse(open(input_fasta), "fasta"))
    for record in records:
        try:
            record.seq = record.seq[nail_start_end[record.id][0]:nail_start_end[record.id][1]]
        except IndexError:
            print("No start or end position found for record:", record.id, "skipping")
            continue
    SeqIO.write(records, output, "fasta")




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("--nail",required=True, help="Number of Ns allowed in the beginning", default=0)
    parser.add_argument("-o", "--output", default="data/nail_start_end_trimmed.fasta", help="Output FASTA file")
    args = parser.parse_args()
    main(args.input, args.nail, args.output)
