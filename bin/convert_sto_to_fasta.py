from Bio import SeqIO
import argparse

def main (input_files):

    for files in input_files:
            
        records = SeqIO.parse(files, "stockholm")
        count = SeqIO.write(records, f"{files.replace(".sto","")}.fasta", "fasta")
        print(f"Converted %i records in file {files}" % count)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Convert Stockholm file into Fasta format")
    parser.add_argument("-i", "--input_files", nargs="+", required=True, help="Input stockholm file(s) to be converted into fasta")
    args = parser.parse_args()

    main(args.input_files)
