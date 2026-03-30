from Bio import SeqIO
import argparse


def main(file, output, frame, unduplicate):

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
    parser.add_argument("-o", "--output", help="Output FASTA file")
    parser.add_argument("--frame", default=1, type=int, help="Frame to extract (default: 1)")
    parser.add_argument("--unduplicate", action="store_true", help="Remove duplicate sequences")
    args = parser.parse_args()
    main(args.input, args.output, args.frame, args.unduplicate)