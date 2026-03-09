from Bio import SeqIO
import argparse

def main (input_files):
    fasta_files = []

    for files in input_files:
        fasta_file = f"{files.replace(".sto","")}.fasta"
        fasta_files.append(fasta_file)
            
        records = SeqIO.parse(files, "stockholm")
        count = SeqIO.write(records, fasta_file, "fasta")

        print(f"Converted %i records in file {files}" % count)
        

    for files in fasta_files:
        # replace the headers to include the SSL number and remove all '-' from the sequence
        ssl_number = files.split("/")[-1].split(".")[0]
        records = list(SeqIO.parse(files, "fasta"))
        for record in records:
            if not record.id.endswith(f"{ssl_number}"):
                record.id = record.id.split(" ")[0]
                record.id = f"{record.id}|{ssl_number}"
            record.seq = record.seq.replace("-","")
            record.description = ""
        SeqIO.write(records, files, "fasta")
          

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Convert Stockholm file into Fasta format")
    parser.add_argument("-i", "--input_files", nargs="+", required=True, help="Input stockholm file(s) to be converted into fasta")
    args = parser.parse_args()

    main(args.input_files)
