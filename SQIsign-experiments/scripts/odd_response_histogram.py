#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0

import argparse
import csv
import os
import statistics
from pathlib import Path


def load_rows(path):
    with path.open(newline="") as fh:
        csv_lines = []
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith("#") or stripped.startswith("Random seed:"):
                continue
            csv_lines.append(line)

    if not csv_lines:
        raise ValueError(f"{path} does not contain CSV rows")

    rows = list(csv.DictReader(csv_lines))
    if not rows:
        raise ValueError(f"{path} does not contain data rows")

    return rows


def select_field(fieldnames, candidates):
    for candidate in candidates:
        if candidate in fieldnames:
            return candidate
    raise ValueError(f"missing one of these columns: {', '.join(candidates)}")


def summarize(values):
    return {
        "min": min(values),
        "max": max(values),
        "mean": statistics.fmean(values),
        "median": statistics.median(values),
    }


def print_summary(label, values, precision=2):
    summary = summarize(values)
    print(
        f"{label}: min={summary['min']} max={summary['max']} "
        f"mean={summary['mean']:.{precision}f} median={summary['median']:.{precision}f}"
    )


def compute_ratios(rows, fieldnames, odd_field, total_field):
    if "odd_responses/total_responses" in fieldnames:
        return [float(row["odd_responses/total_responses"]) for row in rows]

    if "odd_ratio" in fieldnames:
        return [float(row["odd_ratio"]) for row in rows]

    if "odd_percent" in fieldnames:
        return [float(row["odd_percent"]) / 100.0 for row in rows]

    ratios = []
    for row in rows:
        total = int(row[total_field])
        ratios.append(0.0 if total == 0 else int(row[odd_field]) / total)
    return ratios


def binned_counts(values, bins):
    counts = [0] * bins
    for value in values:
        index = int(value * bins)
        if index < 0:
            index = 0
        if index >= bins:
            index = bins - 1
        counts[index] += 1
    return counts


def print_histogram(values, bins, width):
    counts = binned_counts(values, bins)
    max_count = max(counts)

    print("ratio_low,ratio_high,count,histogram")
    for index, count in enumerate(counts):
        if count == 0:
            continue
        low = index / bins
        high = (index + 1) / bins
        bar_len = max(1, round(width * count / max_count))
        print(f"{low:.4f},{high:.4f},{count},{'#' * bar_len}")


def save_histogram(values, path, bins):
    mpl_config = Path("experiments/.matplotlib")
    mpl_config.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config))

    try:
        import matplotlib
    except ModuleNotFoundError as exc:
        if exc.name != "matplotlib":
            raise
        raise SystemExit(
            "matplotlib is not installed for this Python interpreter.\n"
            "Install it in the same environment you use to run this script, for example:\n"
            "  python3 -m pip install matplotlib\n"
            "Or run without --plot to print the text histogram only."
        ) from exc

    matplotlib.use("Agg")
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "matplotlib.pyplot could not be imported. Reinstall matplotlib in this Python environment:\n"
            "  python3 -m pip install --force-reinstall matplotlib"
        ) from exc

    fig, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)
    bin_edges = [i / bins for i in range(bins + 1)]

    ax.hist(values, bins=bin_edges, rwidth=0.9, color="#2f6f9f", edgecolor="#1f2933")
    ax.set_title("Odd response ratio per random commitment/challenge pair")
    ax.set_xlabel("odd_responses / total_responses")
    ax.set_ylabel("frequency")
    ax.set_xlim(0.0, 1.0)
    ax.grid(axis="y", alpha=0.25)

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Print and optionally plot a histogram of odd-response ratios for random commitment/challenge pairs."
    )
    parser.add_argument(
        "csv_file",
        nargs="?",
        default="experiments/lvl1_odd_response_counts.csv",
        help="Counts CSV, or raw odd_response_experiment output. Default: %(default)s",
    )
    parser.add_argument("--width", type=int, default=50, help="Maximum text bar width. Default: %(default)s")
    parser.add_argument("--bins", type=int, default=40, help="Number of ratio bins on [0, 1]. Default: %(default)s")
    parser.add_argument(
        "--plot",
        nargs="?",
        const="",
        default=None,
        help="Save a matplotlib histogram PNG. Optional path defaults to <csv stem>_histogram.png",
    )
    args = parser.parse_args()

    path = Path(args.csv_file)
    rows = load_rows(path)
    fieldnames = rows[0].keys()
    odd_field = select_field(fieldnames, ("odd_responses", "odd_degree"))
    total_field = select_field(fieldnames, ("total_responses", "responses"))

    odd_values = [int(row[odd_field]) for row in rows]
    total_values = [int(row[total_field]) for row in rows]
    ratio_values = compute_ratios(rows, fieldnames, odd_field, total_field)

    print(f"samples: {len(rows)}")
    print_summary("total_responses", total_values)
    print_summary("odd_responses", odd_values)
    print_summary("odd_responses/total_responses", ratio_values, precision=4)
    print_histogram(ratio_values, max(1, args.bins), max(1, args.width))

    if args.plot is not None:
        if args.plot:
            plot_path = Path(args.plot)
        else:
            plot_path = path.with_name(f"{path.stem}_histogram.png")
        save_histogram(ratio_values, plot_path, max(1, args.bins))
        print(f"plot: {plot_path}")


if __name__ == "__main__":
    main()
