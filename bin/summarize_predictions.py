# accession feature score
import argparse
import csv
import os
import glob
from collections import defaultdict
import statistics
import logging
import sys
from pathlib import Path


def setup_logging(verbose: bool, output_dir: Path):
    """Configure logging to console and file."""
    log_level = logging.DEBUG if verbose else logging.INFO
    log_file = output_dir / "jackhmmer_batch.log"  

    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file, mode='w') if verbose else logging.NullHandler(),
            logging.StreamHandler(sys.stdout)
        ]
    )
    if verbose:
        logging.info("Logging initialized. Log file: %s", log_file)

def parse_algpred(algpred_file):
# Subject,ML Score,MERCI Score,BLAST Score,Hybrid Score,Prediction
# extract the hybrid score, each contained in a 2d array to contain: accession, feature, score
    results = []
    with open(algpred_file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0] == 'Subject':
                continue
            results.append({
                'accession': row[0],
                'feature': 'algpred2_hybrid_score',
                'score': float(row[4])
            })
    return results

def parse_bcell(bcell_file):
# Accession,Residue,BepiPred-3.0 score,BepiPred-3.0 linear epitope score
    data = defaultdict(list)
    results = []
    with open(bcell_file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0] == 'Accession':
                continue
            data[row[0]].append(float(row[2]))

    for acc in data:
        # get the max, min, mean, and median
        results.append({
            'accession': acc,
            'feature': 'bepipred3_max_epitope_score',
            'score': max(data[acc])
        })

        results.append({
            'accession': acc,
            'feature': 'bepipred3_min_epitope_score',
            'score': min(data[acc])
        })

        results.append({
            'accession': acc,
            'feature': 'bepipred3_mean_epitope_score',
            'score': statistics.mean(data[acc])
        })

        results.append({
            'accession': acc,
            'feature': 'bepipred3_median_epitope_score',
            'score': statistics.median(data[acc])
        })
    return results

def parse_ifnepitope(ifnepitope_file):
    # Seq_ID,Pattern_ID,Sequence,ML_Score,BLAST_Score,Total_Score,Prediction
    data = defaultdict(list)
    results = []
    with open(ifnepitope_file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0] == 'Seq_ID':
                continue
            data[row[0]].append(float(row[5]))

    for acc in data:
        # get the max, min, mean, and median
        results.append({
            'accession': acc,
            'feature': 'ifnepitope2_max_score',
            'score': max(data[acc])
        })

        results.append({
            'accession': acc,
            'feature': 'ifnepitope2_min_score',
            'score': min(data[acc])
        })

        results.append({
            'accession': acc,
            'feature': 'ifnepitope2_mean_score',
            'score': statistics.mean(data[acc])
        })

        results.append({
            'accession': acc,
            'feature': 'ifnepitope2_median_score',
            'score': statistics.median(data[acc])
        })

    return results

def parse_mhc_dir(directory):
    """
    Parse MHC epitope prediction files in the specified directory.
    Expected files are JSONs with MHC I/II prediction results.
    Filters peptides based on %Rank thresholds for SBs and WBs.
    Returns a list of dictionaries with feature values.
    Args:
        directory (str): Path to the directory containing MHC prediction files.
    Returns:
        List of dictionaries, where each dictionary represents an MHC feature.
    Each dictionary contains:
        - "feature": "mhci" or "mhcii" based on the directory name
        - "subfeature": "score", "percentile", "peptide_length", or "num_peptides"
        - "value": numerical score, length, or count value
    """
    prefix = "mhcii" if "mhcii" in directory else "mhci"
    logging.info(f"Parsing MHC dir {directory} with prefix {prefix}")
    results = []
    try:
        files = [f for f in os.listdir(directory) if not f.endswith(".xls")]
        logging.debug(f"Found {len(files)} out files in MHC dir")
    except Exception as e:
        logging.error(f"Failed listing directory {directory}: {e}")
        return results

    for file in files:
        path = os.path.join(directory, file)
        try:
            with open(path) as f:
                data = f.readlines()
                num_peptides = defaultdict(int)
                scores = defaultdict(list)
                percentiles = defaultdict(list)
                num_sb = defaultdict(int)
                num_wb = defaultdict(int)

                for i, line in enumerate(data):
                    if line.startswith("#") or line.startswith("-") or not line.strip() or "Pos" in line:
                        continue
                    # skip all lines without <=
                    if "<=" not in line:
                        continue

                    parts = line.split()
                    if len(parts) < 10:
                        logging.debug(f"Skipping malformed line {i} in {file}: {line.strip()}")
                        continue
                    try:
                        id = parts[10 if prefix == "mhci" else 7]
                        accession = f'{id.split("_")[0]}'
                        score = float(parts[11 if prefix == "mhci" else 8])
                        percentile = float(parts[12 if prefix == "mhci" else 9])
                        binding_strength = parts[14 if prefix == "mhci" else 12] if len(parts) > (13 if prefix == "mhci" else 11) else "NA"

                        num_peptides[accession] += 1
                        # Filter for Strong Binders only
                        if "SB" in binding_strength:
                            num_sb[accession] += 1
                        elif "WB" in binding_strength:
                            num_wb[accession] += 1

                        scores[accession].append(score)
                        percentiles[accession].append(percentile)
                    except ValueError as e:
                        logging.debug(f"Skipping line {i} in {file} due to conversion error: {e}")

                for accession in num_peptides:
                    avg_score = float(statistics.mean(scores[accession])) if scores[accession] else 0
                    avg_percentile = float(statistics.mean(percentiles[accession])) if percentiles[accession] else 0
                    results.append({"accession": accession, "feature": f"{prefix}_score", "score": avg_score})
                    results.append({"accession": accession, "feature": f"{prefix}_percentile", "score": avg_percentile})
                    results.append({"accession": accession, "feature": f"{prefix}_num_peptides", "score": num_peptides[accession]})
                    results.append({"accession": accession, "feature": f"{prefix}_num_strong_binders", "score": num_sb[accession]})
                    results.append({"accession": accession, "feature": f"{prefix}_num_weak_binders", "score": num_wb[accession]})

                
        except Exception as e:
            logging.error(f"Failed parsing MHC file {file}: {e}")
    logging.info(f"Completed MHC parsing with {len(results)} results")
    return results


def main():
    parser = argparse.ArgumentParser(description='Summarize predictions from a file')
    parser.add_argument('-i', '--input_directory', type=str, required=True, help='Input directory containing prediction files')
    parser.add_argument('-o', '--output_file', type=str, default='prediction_summary.csv', help='Output csv name (default: prediction_summary.csv)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose logging')
    args = parser.parse_args()

    setup_logging(args.verbose, Path(args.output_file).parent)

    input_directory = args.input_directory
    output_file = args.output_file

    # algpred_file = (input_directory, 'algpred', '*_algpred.csv')
    # bcell_file = (input_directory, 'bcell', '*', 'raw_output.csv')
    # ifnepitope_file = (input_directory, 'ifnepitope2', '*_ifnepitope2.csv')
    # mhci_file = os.path.join(input_directory, 'mhci','*_mhci_out')
    # mhcii_file = os.path.join(input_directory, 'mhcii','*_mhcii')

    # get the files according to the pattern each should only have one
    algpred_files = glob.glob(os.path.join(input_directory, 'algpred', '*_algpred.csv'))
    if not algpred_files:
        raise FileNotFoundError("No algpred files found in input directory")
    algpred_file = algpred_files[0]
    

    bcell_files = glob.glob(os.path.join(input_directory,'bcell','*','raw_output.csv'))
    if not bcell_files:
        raise FileNotFoundError("No bcell files found in input directory")
    bcell_file = bcell_files[0]

    ifnepitope_files = glob.glob(os.path.join(input_directory, 'ifnepitope2', '*_ifnepitope2.csv'))
    if not ifnepitope_files:
        raise FileNotFoundError("No ifnepitope files found in input directory")
    ifnepitope_file = ifnepitope_files[0]

    mhci_dir = os.path.join(input_directory, 'mhci')
    if not os.path.exists(mhci_dir):
        raise FileNotFoundError("No mhci files found in input directory")

    mhcii_dir = os.path.join(input_directory, 'mhcii')
    if not os.path.exists(mhcii_dir):
        raise FileNotFoundError("No mhcii files found in input directory")


    algpred_summary = parse_algpred(algpred_file)
    bcell_summary = parse_bcell(bcell_file)
    ifnepitope_summary = parse_ifnepitope(ifnepitope_file)
    mhci_summary = parse_mhc_dir(mhci_dir)
    mhcii_summary = parse_mhc_dir(mhcii_dir)

    all_results = algpred_summary + bcell_summary + ifnepitope_summary + mhci_summary + mhcii_summary

    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["accession", "feature", "score"])
        for result in all_results:
            writer.writerow([result['accession'], result['feature'], result['score']])

    logging.info(f"Results written to {output_file}")


if __name__ == '__main__':
    main()