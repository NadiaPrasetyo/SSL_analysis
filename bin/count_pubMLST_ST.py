import argparse
import csv
import sys
from collections import defaultdict
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np


def main(input_file, output_csv, output_plot=None, year_column=-2, st_column=-1):
    # >3532|Saitama9|Staphylococcus_aureus|Japan|2012|1558
    records = SeqIO.parse(input_file, "fasta")

    # Track counts as {year: {ST: count}}
    year_st_counts = defaultdict(lambda: defaultdict(int))

    for record in records:
        parts = record.id.split("|")
        if len(parts) < 6:
            continue
        year = parts[year_column]
        st = parts[st_column]
        year_st_counts[year][st] += 1

    # Write CSV: year, ST, count
    with open(output_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["year", "ST", "count"])
        for year in sorted(year_st_counts):
            for st, count in sorted(year_st_counts[year].items(), key=lambda x: -x[1]):
                writer.writerow([year, st, count])
    print(f"Wrote {output_csv}")

    # --- Plot ---
    # Find top 5 STs per year (by count in that year)
    years = sorted(year_st_counts.keys())

    # Collect all top-5 STs across all years so we can assign consistent colors
    top_sts_per_year = {}
    all_top_sts = set()
    for year in years:
        top5 = sorted(year_st_counts[year], key=lambda s: -year_st_counts[year][s])[:5]
        top_sts_per_year[year] = top5
        all_top_sts.update(top5)

    # Assign a color to each unique ST
    all_top_sts = sorted(all_top_sts)
    cmap = cm.get_cmap("tab20", len(all_top_sts))
    st_colors = {st: cmap(i) for i, st in enumerate(all_top_sts)}

    fig, ax = plt.subplots(figsize=(12, 6))

    # For each ST that appears in any year's top 5, plot its trajectory
    # Only connect years where it is in the top 5
    for st in all_top_sts:
        x_vals, y_vals = [], []
        for year in years:
            if st in top_sts_per_year[year]:
                x_vals.append(year)
                y_vals.append(year_st_counts[year][st])
        if x_vals:
            ax.plot(x_vals, y_vals, marker="o", label=f"ST {st}",
                    color=st_colors[st], linewidth=2, markersize=6)

    ax.set_xlabel("Year", fontsize=12)
    ax.set_ylabel("ST Count", fontsize=12)
    ax.set_title("Top 5 ST Counts per Year", fontsize=14)
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=9, title="Sequence Type")
    ax.grid(axis="y", linestyle="--", alpha=0.5)
    plt.xticks(rotation=45)
    plt.tight_layout()

    plot_path = output_plot or output_csv.replace(".csv", "_plot.png")
    plt.savefig(plot_path, dpi=150)
    print(f"Wrote {plot_path}")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", type=str, required=True,
                        help="Input file in FASTA format")
    parser.add_argument("-o", "--output_csv", type=str, required=True,
                        help="Output CSV file")
    parser.add_argument("-p", "--output_plot", type=str, default=None,
                        help="Output plot file (default: <output_csv>_plot.png)")
    parser.add_argument("--year-column", type=int, default=-2,
                        help="Column number for year (default: -2 (last column - 1))")
    parser.add_argument("--st-column", type=int, default=-1,
                        help="Column number for ST (default: -1 (last column))")
    args = parser.parse_args()
    main(args.input_file, args.output_csv, args.output_plot, args.year_column, args.st_column)