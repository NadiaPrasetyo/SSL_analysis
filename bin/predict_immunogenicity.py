"""
epitope_pipeline.py

Unified pipeline for:
1. Antigen prediction tools
2. IEDB epitope prediction tools

Supports running:
- antigen tools only
- IEDB tools only
- both pipelines sequentially
"""

import argparse
import sys
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from collections import Counter
from typing import List, Tuple, Callable

from tools import (
    run_algpred,
    run_ifnepitope2,
    run_mhci,
    run_mhcii,
    run_bcell,
    common
)

# --------------------------------------------------
# Tool Registries
# --------------------------------------------------

ANTIGEN_RUNNERS = {
    "ALGPRED": run_algpred.run,
    "IFNEPITOPE2": run_ifnepitope2.run,
}

IEDB_RUNNERS = {
    "MHCI": run_mhci.run,
    "MHCII": run_mhcii.run,
    "BCell": run_bcell.run,
}

Job = Tuple[str, Callable, tuple]


# --------------------------------------------------
# Logging
# --------------------------------------------------

def setup_logging(verbose: bool, output_dir: Path) -> None:
    """Configure logging to console and optional file."""
    output_dir.mkdir(parents=True, exist_ok=True)

    level = logging.DEBUG if verbose else logging.INFO
    log_file = output_dir / "predict_immunogenicity.log"

    handlers = [logging.StreamHandler(sys.stdout)]

    if verbose:
        handlers.append(logging.FileHandler(log_file, mode="w"))

    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=handlers,
    )

    if verbose:
        logging.info("Logging initialized → %s", log_file)


# --------------------------------------------------
# Job Execution
# --------------------------------------------------

def run_tool(tool_name: str, runner: Callable, *args) -> None:
    """Execute a single tool job."""
    try:
        runner(*args)
        logging.info("✅ %s completed", tool_name)
    except Exception as e:
        logging.exception("❌ %s failed: %s", tool_name, e)


def run_parallel_jobs(jobs: List[Job], threads: int) -> None:
    """Run jobs using thread pool."""
    logging.info("⚙️ Running %d jobs using %d threads", len(jobs), threads)

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(run_tool, name, func, *args)
            for name, func, args in jobs
        ]

        for f in futures:
            f.result()


# --------------------------------------------------
# Antigen Pipeline
# --------------------------------------------------

def is_antigen_job_done(fasta: Path, output_dir: Path) -> bool:
    """Check if antigen job already completed."""
    identifier = fasta.stem
    return any(identifier in f.name for f in output_dir.glob("*"))


def build_antigen_jobs(
    args,
    fasta_files: List[Path],
    tool_paths: dict,
    output_root: Path
) -> List[Job]:

    jobs: List[Job] = []

    for tool_name, runner in ANTIGEN_RUNNERS.items():

        if tool_name not in args.tools:
            continue

        output_dir = output_root / tool_name.lower()
        output_dir.mkdir(parents=True, exist_ok=True)

        for fasta in fasta_files:

            if is_antigen_job_done(fasta, output_dir):
                logging.info("⏭️ Skipping %s for %s", tool_name, fasta.name)
                continue

            jobs.append((
                tool_name,
                runner,
                (tool_paths[tool_name], fasta, output_dir, args.batch_size)
            ))

    return jobs


# --------------------------------------------------
# IEDB Pipeline
# --------------------------------------------------

def is_iedb_job_done(tool: str, fasta: Path, output_root: Path) -> bool:
    """Check if IEDB job already completed."""
    subdir = output_root / tool.lower()

    if not subdir.exists():
        return False

    return any(fasta.stem in f.name for f in subdir.glob("*"))


def build_iedb_jobs(
    fasta_files: List[Path],
    tool_map: dict,
    output_root: Path
) -> List[Job]:

    jobs: List[Job] = []

    for tool_name, tool_path in tool_map.items():

        if tool_name not in IEDB_RUNNERS:
            continue

        runner = IEDB_RUNNERS[tool_name]

        for fasta in fasta_files:

            if is_iedb_job_done(tool_name, fasta, output_root):
                logging.info("⏭️ Skipping %s for %s", tool_name, fasta.name)
                continue

            jobs.append((
                tool_name,
                runner,
                (fasta, tool_path, output_root)
            ))

    return jobs


# --------------------------------------------------
# Argument Parsing
# --------------------------------------------------

def parse_args():

    parser = argparse.ArgumentParser(
        description="Unified Antigen + IEDB Epitope Prediction Pipeline"
    )

    parser.add_argument(
        "-i",
        "--input-fasta",
        nargs="+",
        required=True,
        help="Input FASTA file(s)"
    )

    parser.add_argument(
        "--tool-root",
        required=True,
        help="Root directory containing prediction tools"
    )

    parser.add_argument(
        "--pipeline",
        choices=["analysis", "iedb", "all"],
        default="all",
        help="Pipeline to run"
    )

    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--batch-size", type=int, default=10000)

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("immunogenicity"),
        help="Output directory"
    )

    parser.add_argument("--verbose", action="store_true")

    return parser.parse_args()


# --------------------------------------------------
# Main
# --------------------------------------------------

def main():

    args = parse_args()
    setup_logging(args.verbose, args.output_dir)

    fasta_files = [Path(f) for f in args.input_fasta]

    for f in fasta_files:
        if not f.exists():
            logging.error("❌ FASTA file not found: %s", f)
            sys.exit(1)

    tool_root = Path(args.tool_root).resolve()
    output_root = args.output_dir
    output_root.mkdir(parents=True, exist_ok=True)

    all_jobs: List[Job] = []

    # --------------------------------------------------
    # Antigen tools
    # --------------------------------------------------

    if args.pipeline in ("analysis", "all"):

        logging.info("🧪 Preparing Antigen Prediction Jobs")

        tool_paths = common.check_antigen_tools(
            list(ANTIGEN_RUNNERS.keys()),
            tool_root
        )

        all_jobs.extend(
            build_antigen_jobs(
                args,
                fasta_files,
                tool_paths,
                output_root
            )
        )

    # --------------------------------------------------
    # IEDB tools
    # --------------------------------------------------

    if args.pipeline in ("iedb", "all"):

        logging.info("🧬 Preparing IEDB Prediction Jobs")

        tool_map = common.check_iedb_tool(tool_root)

        all_jobs.extend(
            build_iedb_jobs(
                fasta_files,
                tool_map,
                output_root
            )
        )

    # --------------------------------------------------
    # Run jobs
    # --------------------------------------------------

    if not all_jobs:
        logging.info("✅ No jobs to run")
        sys.exit(0)

    job_counter = Counter(job[0] for job in all_jobs)

    logging.info("\n📋 Job Summary:")
    for tool, count in job_counter.items():
        logging.info("  %s: %d", tool, count)

    run_parallel_jobs(all_jobs, args.threads)

    logging.info("\n✅ Pipeline complete")


if __name__ == "__main__":
    main()