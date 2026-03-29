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
    1.  Strip any FASTA-style description after the first whitespace.
        "ACL82857.1 carbamate kinase partial [...]"  →  "ACL82857.1"

    2a. "ACL82857|1"  →  "ACL82857.1"
        A pipe whose RIGHT side is a bare integer AND whose left side
        contains no underscore (plain accession + version).  The pipe
        is treated as a mis-encoded dot.

    2b. "YP|499715"  →  "YP_499715"
        A pipe whose LEFT side is a short all-letter prefix (2-4 chars)
        AND whose right side is all-digits with no trailing version.
        This is an NCBI-style accession where '_' was encoded as '|'
        by the MHC prediction tool.  We restore the underscore.
        The version suffix (.1) is intentionally left absent here so that
        resolve_accession can match it via the version-stripped variant.

    Rule 2b is tested BEFORE 2a so that "YP|499715" is caught first.

    3.  Everything else containing '|' (e.g. "SSL3|CC1", "SSL10|A0A0H3K6Z9",
        "sp|P64362.1|HIS6_STAAN") is left unchanged.
    """
    # Step 1 – drop trailing description
    acc = raw.strip().split()[0] if raw.strip() else raw.strip()

    # Step 2b – "YP|499715" → "YP_499715"
    # Left side: 2-4 uppercase letters; right side: digits only (no dot)
    acc = re.sub(r'^([A-Z]{2,4})\|(\d+)$', r'\1_\2', acc)

    # Step 2a – "ACL82857|1" → "ACL82857.1"  (trailing |<integer> = version)
    acc = re.sub(r'\|(\d+)$', r'.\1', acc)

    return acc


def build_accession_map(canonical_accessions: list[str]):
    """
    Build a lookup dict from every plausible variant of each canonical
    accession to the canonical form itself.

    For every canonical we register all the mangled forms that real
    prediction tools have been observed to produce:

    Canonical "YP_499715.1":
      YP_499715.1          – identity
      yp_499715.1          – lower-cased
      YP_499715|1          – dot→pipe  (some tools)
      YP|499715.1          – underscore→pipe
      YP|499715            – underscore→pipe + version stripped  ← MHC format
      YP_499715            – version stripped (no .1)
      YP.499715.1          – underscore→dot  (occasionally seen)
      YP.499715            – underscore→dot + version stripped

    Canonical "ACL82857.1":
      ACL82857.1  ACL82857|1  acl82857.1  ACL82857  …

    Canonical "SSL3|CC1":
      SSL3|CC1  ssl3|cc1  SSL3.CC1  …  (pipe IDs kept as-is)

    Canonical "sp|P64362.1|HIS6_STAAN":
      kept as-is plus lower-cased; we do not mangle multi-pipe UniProt IDs
      further to avoid collisions.
    """
    mapping: dict[str, str] = {}

    for canon in canonical_accessions:
        variants: set[str] = set()
        variants.add(canon)
        variants.add(canon.lower())

        # Only generate separator-swap variants for simple accessions
        # (at most one underscore separator, at most one dot version suffix).
        # Multi-pipe UniProt IDs like "sp|P64362.1|HIS6_STAAN" are kept verbatim.
        pipe_count = canon.count('|')
        has_underscore = '_' in canon

        if pipe_count <= 1:
            # dot ↔ pipe for the version suffix
            variants.add(canon.replace('.', '|'))   # YP_499715|1
            variants.add(canon.replace('|', '.'))   # SSL3.CC1
            # pipe → underscore (MHC tools sometimes encode SSL10|A0A0H3K6Z9 as SSL10_A0A0H3K6Z9)
            variants.add(canon.replace('|', '_'))   # SSL10_A0A0H3K6Z9

        if has_underscore and pipe_count == 0:
            # NCBI-style: swap underscore for pipe
            us_to_pipe = canon.replace('_', '|')    # YP|499715.1
            variants.add(us_to_pipe)

            # MHC tool drops the .version entirely after the swap
            us_pipe_no_ver = re.sub(r'\.\d+$', '', us_to_pipe)  # YP|499715
            variants.add(us_pipe_no_ver)

            # Also dot-separated variants (seen in auto-registered accessions)
            us_to_dot = canon.replace('_', '.')     # YP.499715.1
            variants.add(us_to_dot)
            variants.add(re.sub(r'\.\d+$', '', us_to_dot))      # YP.499715

        # Version-stripped form (drop trailing .N)
        no_version = re.sub(r'\.\d+$', '', canon)
        if no_version != canon:
            variants.add(no_version)
            variants.add(no_version.lower())

        for v in variants:
            if v and (v not in mapping or len(v) >= len(mapping.get(v, ''))):
                mapping[v] = canon

    # Build a secondary index keyed by separator-flattened canonical.
    # This lets resolve_accession handle truncated all-underscore MHC ids.
    # We index by the flattened *canonical* (not variants) so that the
    # startswith probe in resolve_accession matches the longest real prefix.
    flat_index: dict[str, str] = {}
    for canon in canonical_accessions:
        flat_canon = _flatten(canon)
        # Only add if not already claimed by a longer canonical
        if flat_canon not in flat_index or len(flat_canon) >= len(_flatten(flat_index[flat_canon])):
            flat_index[flat_canon] = canon

    return mapping, flat_index


def _flatten(s: str) -> str:
    """Replace all separator characters (|, ., -) with _ for fuzzy comparison."""
    return re.sub(r'[|.\-]', '_', s)


def resolve_accession(raw: str, acc_map: dict[str, str],
                      flat_index: dict[str, str]) -> str:
    """
    Return the canonical accession for *raw*, using acc_map.

    MHC prediction tools mangle accession IDs in two ways that must both
    be handled:
      - All separators (|  .  -) are replaced with underscores
      - The result is truncated to a fixed column width

    Lookup order:
    1. Exact match on the normalised string
    2. Exact match on the raw first token (catches pre-normalised keys)
    3. Case-insensitive match on normalised string
    4. Separator-variant prefix match via acc_map variants
       (catches YP_499715 → YP_499715.1, SSL10_A0A0H3K6Z9 → SSL10|A0A0H3K6Z9)
    5. Separator-flattened prefix match via flat_index
       (catches truncated all-underscore forms like:
        SSL10_A0A0H3K6Z  → SSL10|A0A0H3K6Z9
        P0C0I5_Exotoxin  → P0C0I5|ExotoxinC
        sp_P64362_1_HIS  → sp|P64362.1|HIS6_STAAN)
    6. Unknown – auto-register for cross-tool consistency.
    """
    norm    = normalize_accession(raw)
    raw_tok = raw.strip().split()[0]

    # 1. Exact on normalised
    if norm in acc_map:
        return acc_map[norm]

    # 2. Exact on raw token
    if raw_tok in acc_map:
        return acc_map[raw_tok]

    # 3. Case-insensitive
    if norm.lower() in acc_map:
        return acc_map[norm.lower()]

    # 4. Longest-prefix match across all registered separator variants
    best, best_len = None, 0
    for variant, canon in acc_map.items():
        if len(variant) >= 6 and norm.startswith(variant) and len(variant) > best_len:
            best, best_len = canon, len(variant)
    if best:
        return best

    # 5. Separator-flattened prefix match.
    #    The MHC tool truncates the id, so the probe is SHORTER than the
    #    canonical.  We check: does the flat canonical START WITH the probe?
    #    Longest probe match wins (avoids short spurious hits).
    flat_norm = _flatten(norm)
    flat_tok  = _flatten(raw_tok)
    best, best_len = None, 0
    for flat_canon, canon in flat_index.items():
        for probe in (flat_norm, flat_tok):
            if len(probe) >= 6 and flat_canon.startswith(probe) and len(probe) > best_len:
                best, best_len = canon, len(probe)
    if best:
        logging.debug(
            "Resolved %r via flattened-prefix match -> %r", raw, best
        )
        return best

    # 6. Never-seen-before – register for consistency
    logging.debug(
        "New accession %r (normalised: %r) – adding to map automatically.", raw, norm
    )
    new_variants = build_accession_map([norm])
    acc_map.update(new_variants)
    flat_index.update({_flatten(v): norm for v in new_variants})
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

def parse_algpred(algpred_file, acc_map, flat_index):
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
            acc = resolve_accession(row[0], acc_map, flat_index)
            results.append({
                'accession': acc,
                'feature': 'algpred2_hybrid_score',
                'score': score
            })
    return results


def parse_bcell(bcell_file, acc_map, flat_index):
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
            acc = resolve_accession(row[0], acc_map, flat_index)
            data[acc].append(score)

    for acc, scores in data.items():
        results.append({'accession': acc, 'feature': 'bepipred3_max_epitope_score',    'score': max(scores)})
        results.append({'accession': acc, 'feature': 'bepipred3_min_epitope_score',    'score': min(scores)})
        results.append({'accession': acc, 'feature': 'bepipred3_mean_epitope_score',   'score': statistics.mean(scores)})
        results.append({'accession': acc, 'feature': 'bepipred3_median_epitope_score', 'score': statistics.median(scores)})
    return results


def parse_ifnepitope(ifnepitope_file, acc_map, flat_index):
    # Seq_ID,Pattern_ID,Sequence,ML_Score,BLAST_Score,Total_Score,Prediction
    data = defaultdict(list)
    results = []
    with open(ifnepitope_file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0] == 'Seq_ID':
                continue
            acc = resolve_accession(row[0], acc_map, flat_index)
            try:
                data[acc].append(float(row[5]))
            except (ValueError, IndexError):
                logging.debug("Skipping ifnepitope row (could not parse score): %s", row)
                continue

    for acc, scores in data.items():
        results.append({'accession': acc, 'feature': 'ifnepitope2_max_score',    'score': max(scores)})
        results.append({'accession': acc, 'feature': 'ifnepitope2_min_score',    'score': min(scores)})
        results.append({'accession': acc, 'feature': 'ifnepitope2_mean_score',   'score': statistics.mean(scores)})
        results.append({'accession': acc, 'feature': 'ifnepitope2_median_score', 'score': statistics.median(scores)})
    return results


def parse_mhc_dir(directory, acc_map, flat_index):
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
                        # MHC files encode accessions as underscore-joined tokens,
                        # e.g. "SSL10_A0A0H3K6Z9" or "YP_499715" (underscore used
                        # instead of pipe/dot).  Pass raw_id directly to
                        # resolve_accession; normalize_accession and the variant
                        # map handle all the separator permutations correctly.
                        # We do NOT split and rejoin here because that would
                        # truncate IDs with multiple underscores (e.g. A0A0H3K6Z9).
                        accession = resolve_accession(raw_id, acc_map, flat_index)

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

    acc_map, flat_index = build_accession_map(canonicals)

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
    algpred_summary    = parse_algpred(algpred_file,      acc_map, flat_index)
    bcell_summary      = parse_bcell(bcell_file,           acc_map, flat_index)
    ifnepitope_summary = parse_ifnepitope(ifnepitope_file, acc_map, flat_index)
    mhci_summary       = parse_mhc_dir(mhci_dir,          acc_map, flat_index)
    mhcii_summary      = parse_mhc_dir(mhcii_dir,         acc_map, flat_index)

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