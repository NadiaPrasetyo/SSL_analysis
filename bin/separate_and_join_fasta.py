from Bio import SeqIO
import argparse
from pathlib import Path

def join_fasta(input_dir):
    """
    Concatenate all FASTA files in a directory into one list of records.

    Parameters
    ----------
    input_dir : Path
        Directory containing FASTA files to concatenate.

    Returns
    -------
    list[SeqRecord]
        A list of SeqRecord objects containing the concatenated FASTA records.
    """
    print(f"Joining FASTA files from {input_dir}")
    records = []

    for file in input_dir.iterdir():
        records.extend(SeqIO.parse(file, "fasta"))

    return records

def separate_fasta(records, output_dir):
    print(f"Separating FASTA records into {output_dir}")
    num_SSL3 = 0
    num_SSL7 = 0
    num_SSL11 = 0
    for record in records:
        # >3882|2505034|ERR410047|Unknown|2016|ST-39|Q2G0X4|SSL11 ali_from=16582 ali_to=17274 evalue=1.9e-82
        record.id = record.id.split(" ")[0]
        print(record.id)
        if len(record.id.split("|")) >= 8:
            record.id = record.id.split("|")[2] + "_" + record.id.split("|")[3] + "_" + record.id.split("|")[4] + "_" + (record.id.split("|")[5].replace("-","")) + "_" + record.id.split("|")[7]
        record.description = ""
        if "SSL3" in record.id:
            num_SSL3 += 1
            with open(output_dir / "SSL3.fasta", "a") as f:
                SeqIO.write(record, f, "fasta")
        elif "SSL7" in record.id:
            num_SSL7 += 1
            with open(output_dir / "SSL7.fasta", "a") as f:
                SeqIO.write(record, f, "fasta")
        elif "SSL11" in record.id:
            num_SSL11 += 1
            with open(output_dir / "SSL11.fasta", "a") as f:
                SeqIO.write(record, f, "fasta")

        else:
            with open(output_dir / "not_matched.fasta", "a") as f:
                SeqIO.write(record, f, "fasta")

    print(f"Separated FASTA records into {output_dir}, number of records:\nSSL3: {num_SSL3}\nSSL7: {num_SSL7}\nSSL11: {num_SSL11}")
    return

def filter_year(records, year):
    new_records = []
    for record in records:
        if record.id.split("|")[4] >= year:
            new_records.append(record)
    return new_records

def deduplicate_fasta_records(records):
    seen_seqs = set()
    deduped_records = []
    # check the headers that they have the correct pattern before deduplicating
    for record in records:
        if len(record.id.split("|")) < 8:
            # skip this record
            continue
        if record.seq not in seen_seqs and not seen_seqs.add(record.seq):
            deduped_records.append(record)
    return deduped_records

def remove_unknown_records(records):
    print("Removing records with 'unknown' in the header")
    new_records = []
    for record in records:
        if "unknown" not in record.id.lower():
            new_records.append(record)
    return new_records

def main():
    parser = argparse.ArgumentParser(description="Separate and join FASTA files based on SSL number.")
    parser.add_argument("-i", "--input-dir", required=True, help="Input directory containing FASTA files.")
    parser.add_argument("-o", "--output-dir", default="separated_fasta", help="Output directory for separated FASTA files. (default: separated_fasta)")
    parser.add_argument("--deduplicate", action="store_true", help="Deduplicate the FASTA records.")
    parser.add_argument("--remove-unknown", action="store_true", help="Remove unknown records.")
    parser.add_argument("--filter-year", type=int, help="Filter records by year. (default: None)")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)
    
    records = join_fasta(input_dir)

    if args.remove_unknown:
        print(f"Removing unknown records: starting with {len(records)} records.")
        records = remove_unknown_records(records)
        print(f"Removed unknown records: {len(records)} records.")

    if args.deduplicate:
        print(f"Deduplicating FASTA records: starting with {len(records)} records.")
        records = deduplicate_fasta_records(records)
        print(f"Deduplicated FASTA records: {len(records)} records.")

    if args.filter_year != None:
        print(f"Filtering records by year: starting with {len(records)} records.")
        records = filter_year(records, args.filter_year)
        print(f"Filtered records by year: {len(records)} records.")

    ref_path = "data/SSL_seq.fasta"
    chosen_path = "data/6_strains/aa/SSL_aa.fasta"

    print(f"Adding reference and chosen sequences from {ref_path} and {chosen_path}")

    # Add reference and chosen sequences to the records
    records.extend(SeqIO.parse(ref_path, "fasta"))
    records.extend(SeqIO.parse(chosen_path, "fasta"))

    separate_fasta(records, output_dir)

if __name__ == "__main__":
    main()
