from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import csv

def tsv_file_to_fasta(tsv_file):
    # tsv header:
    # query	target	pident	alnlen	evalue	bits	qcov	tcov	tstart	tend	taln
    seen_ids = set()
    records = []
    with open(tsv_file, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if row[0] == "query":
                continue

            if row[1] in seen_ids:
                continue
                
            seen_ids.add(row[1])
            record = SeqIO.SeqRecord(
                Seq(row[10].replace("-", "")), id=row[1], description=""
            )
            records.append(record)

    print("Number of sequences in tsv file: %i" % len(records)) # print number of records
    return records


def main(file, tsv_file, output, frame, unduplicate):
    if tsv_file is not None:
        records = tsv_file_to_fasta(tsv_file)
    else:
        records = list(SeqIO.parse(file, "fasta"))
    frame_records = []

    for record in records:
        if record.id.endswith(f"frame={frame}"):
            frame_records.append(record)

        if unduplicate:
            # remove duplicate sequences
            seen_seqs = set()
            frame_records_undup = [record for record in frame_records if record.seq not in seen_seqs and not seen_seqs.add(record.seq)]

            with open(output, "w") as f:
                SeqIO.write(frame_records_undup, f, "fasta")
            print("Number of unique sequences in frame %i: %i" % (frame, len(frame_records_undup)))

        else:
            with open(output, "w") as f:
                SeqIO.write(frame_records, f, "fasta")
                
            print("Number of sequences in frame %i: %i" % (frame, len(frame_records)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i ", "--input", help="Input FASTA file")
    parser.add_argument("--input-tsv", help="Input tsv file")
    parser.add_argument("-o", "--output", default="data/extracted_frame_n.fasta",help="Output FASTA file")
    parser.add_argument("--frame", default=1, type=int, help="Frame to extract (default: 1)")
    parser.add_argument("--unduplicate", action="store_true", help="Remove duplicate sequences")
    args = parser.parse_args()

    # check that input or input tsv has been provided
    if args.input is None and args.input_tsv is None:
        raise ValueError("Please provide either an input FASTA file or an input tsv file")

    main(args.input, args.input_tsv, args.output, args.frame, args.unduplicate)