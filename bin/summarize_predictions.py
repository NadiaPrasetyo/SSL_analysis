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
import re


# ---------------------------------------------------------------------------
# Accession normalisation
# ---------------------------------------------------------------------------

def normalize_accession(raw: str) -> str:
    """
    Return a canonical accession string so that the same protein is
    represented identically regardless of which tool produced the label.

    Rules applied in order:
    1. Strip any FASTA-style description after the first whitespace
       e.g. "ACL82857.1 carbamate kinase partial [...]"  →  "ACL82857.1"
    2. Replace pipe-separated version numbers that look like a dot-version
       e.g. "ACL82857|1"  →  "ACL82857.1"
       but leave deliberate pipe IDs like "SSL3|CC1" alone (second part is
       non-numeric or contains letters beyond a simple integer).
    3. Collapse any remaining runs of whitespace.
    """
    # Step 1 – drop trailing description
    acc = raw.strip().split()[0] if raw.strip() else raw.strip()

    # Step 2 – "ACL82857|1" → "ACL82857.1"  (pipe followed by a pure integer)
    acc = re.sub(r'\|(\d+)$', r'.\1', acc)

    return acc


def build_accession_map(canonical_accessions: list[str]):
    """
    Build a lookup dict from every plausible variant of each canonical
    accession to the canonical form itself, so that downstream matching
    can be done with a single dict.get() call.

    Variants generated for each canonical (e.g. "YP_499715.1"):
    - The canonical string itself
    - Lower-cased version
    - '.' replaced by '|'           →  "YP_499715|1"
    - '|' replaced by '.'           →  (handles SSL3|CC1 → SSL3.CC1 lookup)
    - '_' replaced by '|'           →  "YP|499715.1"
    - '_' replaced by '|', '.' → '' →  "YP|499715"   (MHC file format)
    - bare prefix before first separator →  "YP_499715", "ACL82857"
    """
    mapping: dict[str, str] = {}

    for canon in canonical_accessions:
        variants = set()
        variants.add(canon)
        variants.add(canon.lower())
        variants.add(canon.replace('.', '|'))
        variants.add(canon.replace('|', '.'))
        variants.add(canon.replace('_', '|'))
        # e.g. "YP_499715.1" → "YP|499715" (drop .version after underscore swap)
        variants.add(re.sub(r'\.\d+$', '', canon.replace('_', '|')))
        # bare prefix (before the first separator of any kind)
        prefix = re.split(r'[.|_]', canon)[0]
        variants.add(prefix)
        variants.add(prefix.lower())

        for v in variants:
            if v not in mapping or len(v) >= len(mapping.get(v, '')):
                mapping[v] = canon

    return mapping


def resolve_accession(raw: str, acc_map: dict[str, str]) -> str:
    """
    Return the canonical accession for *raw*, using acc_map.

    Lookup order:
    1. Exact match on the normalised string
    2. Case-insensitive match
    3. Prefix match: return the canonical whose variant is a prefix of
       the normalised raw accession (longest prefix wins)
    4. Unknown accession – auto-register all its variants into acc_map
       so that the same protein resolves to the same string across all
       tools even when it wasn't in the original FASTA.
    """
    norm = normalize_accession(raw)

    # 1. Exact
    if norm in acc_map:
        return acc_map[norm]

    # 2. Case-insensitive
    if norm.lower() in acc_map:
        return acc_map[norm.lower()]

    # 3. Prefix: find the canonical whose registered variant is the
    #    longest prefix of norm (guards against short spurious matches)
    best = None
    best_len = 0
    for variant, canon in acc_map.items():
        if len(variant) >= 6 and norm.startswith(variant) and len(variant) > best_len:
            best = canon
            best_len = len(variant)
    if best:
        return best

    # 4. Never-seen-before accession.
    #    Register it now so all subsequent tools resolve it identically.
    logging.debug(
        "New accession %r (normalised: %r) – adding to map automatically.", raw, norm
    )
    # Reuse build_accession_map logic for a single entry
    new_variants = build_accession_map([norm])
    acc_map.update(new_variants)
    return norm


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Per-tool parsers  (each now accepts acc_map and calls resolve_accession)
# ---------------------------------------------------------------------------

def parse_algpred(algpred_file, acc_map):
    # Subject,ML Score,MERCI Score,BLAST Score,Hybrid Score,Prediction
    results = []
    with open(algpred_file, 'r', encoding='utf-8-sig') as f:
        reader = csv.reader(f)
        for row in reader:
            if not row or row[0].strip().lower() in ('subject', ''):
                continue
            try:
                score = float(row[4])
            except (ValueError, IndexError):
                logging.debug("Skipping algpred row (could not parse score): %s", row)
                continue
            acc = resolve_accession(row[0], acc_map)
            results.append({
                'accession': acc,
                'feature': 'algpred2_hybrid_score',
                'score': score
            })
    return results


def parse_bcell(bcell_file, acc_map):
    # Accession,Residue,BepiPred-3.0 score,BepiPred-3.0 linear epitope score
    data = defaultdict(list)
    results = []
    with open(bcell_file, 'r', encoding='utf-8-sig') as f:  # utf-8-sig strips BOM if present
        reader = csv.reader(f)
        for row in reader:
            if not row or row[0].strip().lower() in ('accession', ''):
                continue
            try:
                score = float(row[2])
            except (ValueError, IndexError):
                logging.debug("Skipping bcell row (could not parse score): %s", row)
                continue
            acc = resolve_accession(row[0], acc_map)
            data[acc].append(score)

    for acc, scores in data.items():
        results.append({'accession': acc, 'feature': 'bepipred3_max_epitope_score',    'score': max(scores)})
        results.append({'accession': acc, 'feature': 'bepipred3_min_epitope_score',    'score': min(scores)})
        results.append({'accession': acc, 'feature': 'bepipred3_mean_epitope_score',   'score': statistics.mean(scores)})
        results.append({'accession': acc, 'feature': 'bepipred3_median_epitope_score', 'score': statistics.median(scores)})
    return results


def parse_ifnepitope(ifnepitope_file, acc_map):
    # Seq_ID,Pattern_ID,Sequence,ML_Score,BLAST_Score,Total_Score,Prediction
    data = defaultdict(list)
    results = []
    with open(ifnepitope_file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0] == 'Seq_ID':
                continue
            acc = resolve_accession(row[0], acc_map)
            data[acc].append(float(row[5]))

    for acc, scores in data.items():
        results.append({'accession': acc, 'feature': 'ifnepitope2_max_score',    'score': max(scores)})
        results.append({'accession': acc, 'feature': 'ifnepitope2_min_score',    'score': min(scores)})
        results.append({'accession': acc, 'feature': 'ifnepitope2_mean_score',   'score': statistics.mean(scores)})
        results.append({'accession': acc, 'feature': 'ifnepitope2_median_score', 'score': statistics.median(scores)})
    return results


def parse_mhc_dir(directory, acc_map):
    """
    Parse MHC epitope prediction files in the specified directory.
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
                num_peptides  = defaultdict(int)
                scores        = defaultdict(list)
                percentiles   = defaultdict(list)
                num_sb        = defaultdict(int)
                num_wb        = defaultdict(int)

                for i, line in enumerate(data):
                    if line.startswith("#") or line.startswith("-") or not line.strip() or "Pos" in line:
                        continue
                    if "<=" not in line:
                        continue

                    parts = line.split()
                    if len(parts) < 10:
                        logging.debug(f"Skipping malformed line {i} in {file}: {line.strip()}")
                        continue
                    try:
                        raw_id   = parts[10 if prefix == "mhci" else 7]
                        # MHC files use "GENE_ACCPART1_ACCPART2" style ids;
                        # reconstruct the pipe-separated accession before resolving.
                        id_parts = raw_id.split("_")
                        raw_acc  = f'{id_parts[0]}|{id_parts[1]}' if len(id_parts) >= 2 else raw_id
                        accession = resolve_accession(raw_acc, acc_map)

                        score            = float(parts[11 if prefix == "mhci" else 8])
                        percentile       = float(parts[12 if prefix == "mhci" else 9])
                        binding_strength = parts[14 if prefix == "mhci" else 12] if len(parts) > (13 if prefix == "mhci" else 11) else "NA"

                        num_peptides[accession] += 1
                        if "SB" in binding_strength:
                            num_sb[accession] += 1
                        elif "WB" in binding_strength:
                            num_wb[accession] += 1

                        scores[accession].append(score)
                        percentiles[accession].append(percentile)
                    except (ValueError, IndexError) as e:
                        logging.debug(f"Skipping line {i} in {file} due to error: {e}")

                for accession in num_peptides:
                    avg_score      = float(statistics.mean(scores[accession]))      if scores[accession]      else 0
                    avg_percentile = float(statistics.mean(percentiles[accession])) if percentiles[accession] else 0
                    results.append({"accession": accession, "feature": f"{prefix}_score",               "score": avg_score})
                    results.append({"accession": accession, "feature": f"{prefix}_percentile",          "score": avg_percentile})
                    results.append({"accession": accession, "feature": f"{prefix}_num_peptides",        "score": num_peptides[accession]})
                    results.append({"accession": accession, "feature": f"{prefix}_num_strong_binders",  "score": num_sb[accession]})
                    results.append({"accession": accession, "feature": f"{prefix}_num_weak_binders",    "score": num_wb[accession]})

        except Exception as e:
            logging.error(f"Failed parsing MHC file {file}: {e}")

    logging.info(f"Completed MHC parsing with {len(results)} results")
    return results


# ---------------------------------------------------------------------------
# Canonical accession list  (loaded from the FASTA headers you provide)
# ---------------------------------------------------------------------------

def load_canonical_accessions(fasta_headers: list[str]) -> list[str]:
    """
    Extract and normalise the accession from each FASTA header line.
    Accepts lines with or without the leading '>'.
    """
    canonicals = []
    for line in fasta_headers:
        line = line.strip().lstrip('>')
        if line:
            canonicals.append(normalize_accession(line))
    return canonicals


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description='Summarize predictions from a file')
    parser.add_argument('-i', '--input_directory', type=str, required=True,
                        help='Input directory containing prediction files')
    parser.add_argument('-o', '--output_file', type=str, default='prediction_summary.csv',
                        help='Output csv name (default: prediction_summary.csv)')
    parser.add_argument('-f', '--fasta', type=str, default=None,
                        help='FASTA file used as input (used to build canonical accession list). '
                             'If omitted, accessions are normalised on-the-fly without cross-referencing.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose logging')
    args = parser.parse_args()

    setup_logging(args.verbose, Path(args.output_file).parent)

    input_directory = args.input_directory
    output_file     = args.output_file

    # ------------------------------------------------------------------
    # Build the canonical accession map
    # ------------------------------------------------------------------
    if args.fasta:
        with open(args.fasta) as fh:
            headers = [l for l in fh if l.startswith('>')]
        canonicals = load_canonical_accessions(headers)
        logging.info("Loaded %d canonical accessions from %s", len(canonicals), args.fasta)
    else:
        # We'll collect accessions lazily; normalisation still helps.
        canonicals = []
        logging.info("No FASTA provided – accessions will be normalised but not cross-referenced.")

    acc_map = build_accession_map(canonicals)

    # ------------------------------------------------------------------
    # Locate input files
    # ------------------------------------------------------------------
    algpred_files = glob.glob(os.path.join(input_directory, 'algpred', '*_algpred.csv'))
    if not algpred_files:
        raise FileNotFoundError("No algpred files found in input directory")
    algpred_file = algpred_files[0]

    bcell_files = glob.glob(os.path.join(input_directory, 'bcell', '*', 'raw_output.csv'))
    if not bcell_files:
        raise FileNotFoundError("No bcell files found in input directory")
    bcell_file = bcell_files[0]

    ifnepitope_files = glob.glob(os.path.join(input_directory, 'ifnepitope2', '*_ifnepitope2.csv'))
    if not ifnepitope_files:
        raise FileNotFoundError("No ifnepitope files found in input directory")
    ifnepitope_file = ifnepitope_files[0]

    mhci_dir = os.path.join(input_directory, 'mhci')
    if not os.path.exists(mhci_dir):
        raise FileNotFoundError("No mhci directory found in input directory")

    mhcii_dir = os.path.join(input_directory, 'mhcii')
    if not os.path.exists(mhcii_dir):
        raise FileNotFoundError("No mhcii directory found in input directory")

    # ------------------------------------------------------------------
    # Parse
    # ------------------------------------------------------------------
    algpred_summary  = parse_algpred(algpred_file,    acc_map)
    bcell_summary    = parse_bcell(bcell_file,         acc_map)
    ifnepitope_summary = parse_ifnepitope(ifnepitope_file, acc_map)
    mhci_summary     = parse_mhc_dir(mhci_dir,        acc_map)
    mhcii_summary    = parse_mhc_dir(mhcii_dir,       acc_map)

    all_results = algpred_summary + bcell_summary + ifnepitope_summary + mhci_summary + mhcii_summary

    # ------------------------------------------------------------------
    # Write output
    # ------------------------------------------------------------------
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["accession", "feature", "score"])
        for result in all_results:
            writer.writerow([result['accession'], result['feature'], result['score']])

    logging.info("Results written to %s", output_file)

    # Summarise which accessions appeared
    seen    = {r['accession'] for r in all_results}
    known   = seen & set(canonicals)
    unknown = seen - set(canonicals)
    logging.info("Unique accessions in output : %d", len(seen))
    logging.info("  Matched to FASTA          : %d  ->  %s", len(known), sorted(known))
    if unknown:
        logging.warning(
            "  Not in FASTA (auto-registered, check inputs): %d  ->  %s",
            len(unknown), sorted(unknown)
        )


if __name__ == '__main__':
    main()