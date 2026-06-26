import argparse
import csv
from collections import defaultdict
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def main(input_file, output_csv, output_plot=None, year_column=-2, st_column=-1, plot_year=None, last_n_years=10):
    # >3532|Saitama9|Staphylococcus_aureus|Japan|2012|1558
    records = SeqIO.parse(input_file, "fasta")

    # Track counts as {year: {ST: count}}
    year_st_counts = defaultdict(lambda: defaultdict(int))

    for record in records:
        parts = record.id.split("|")
        if len(parts) < 6:
            continue
        year_raw = parts[year_column].strip()
        st = parts[st_column].strip()

        # Skip missing/empty year or ST
        if not year_raw or not st:
            continue
        try:
            year = int(year_raw)
        except ValueError:
            continue

        year_st_counts[year][st] += 1

    # Write CSV: year, ST, count (sorted by year numerically, then descending count)
    with open(output_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["year", "ST", "count"])
        for year in sorted(year_st_counts):
            for st, count in sorted(year_st_counts[year].items(), key=lambda x: -x[1]):
                writer.writerow([year, st, count])
    print(f"Wrote {output_csv}")

    # --- Plot ---
    all_years = sorted(year_st_counts.keys())  # now sorted numerically

    if plot_year:
        years = [y for y in all_years if y in plot_year]
    elif last_n_years:
        years = all_years[-last_n_years:]
    else:
        years = all_years

    # Find top 5 STs per year
    top_sts_per_year = {}
    all_top_sts = set()
    for year in years:
        top5 = sorted(year_st_counts[year], key=lambda s: -year_st_counts[year][s])[:5]
        top_sts_per_year[year] = top5
        all_top_sts.update(top5)

    # Sort STs and assign consistent colors
    all_top_sts = sorted(all_top_sts, key=lambda s: int(s) if s.isdigit() else float("inf"))
    cmap = cm.get_cmap("tab20", max(len(all_top_sts), 1))
    st_colors = {st: cmap(i) for i, st in enumerate(all_top_sts)}

    fig, ax = plt.subplots(figsize=(14, 6))

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
    ax.grid(axis="y", linestyle="--", alpha=0.5)

    # Set x-axis to only show the plotted years, as integers
    ax.set_xticks(years)
    ax.set_xticklabels([str(y) for y in years], rotation=45, ha="right")
    ax.set_xlim(years[0] - 0.5, years[-1] + 0.5)

    # Legend outside plot, but constrained so it doesn't overflow
    ax.legend(
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        fontsize=8,
        title="Sequence Type",
        title_fontsize=9,
        frameon=True,
        borderpad=0.5,
        ncol=max(1, len(all_top_sts) // 20),  # use 2 columns if many STs
    )

    plt.tight_layout()
    plot_path = output_plot or output_csv.replace(".csv", "_plot.png")
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
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
                        help="Column index for year (default: -2)")
    parser.add_argument("--st-column", type=int, default=-1,
                        help="Column index for ST (default: -1)")
    parser.add_argument("--plot-year", type=int, default=None, nargs="+",
                        help="Specific years to plot (overrides --last-n-years)")
    parser.add_argument("--last-n-years", type=int, default=10,
                        help="Plot only the last N years (default: 10); ignored if --plot-year is set")
    args = parser.parse_args()
    main(args.input_file, args.output_csv, args.output_plot,
         args.year_column, args.st_column, args.plot_year, args.last_n_years)