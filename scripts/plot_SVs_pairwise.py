#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt


def arg_parser():
    parser = argparse.ArgumentParser(
        prog='plot_SVs_pairwise.py',
        description='Plot SV length distributions (insertions/deletions) from BED files.'
    )

    parser.add_argument(
        "-b", "--bed_dir",
        required=True,
        help="Directory containing BED files with columns: sample1 start1 end1 sample2 start2 end2 type"
    )

    parser.add_argument(
        "-o", "--output",
        default="sv_length_violinplot.png",
        help="Output filename for the violin plot (default: sv_length_violinplot.png)"
    )

    return parser.parse_args()


def read_bed_file(filepath):
    """
    Reads a bed-like file with columns:
    sample1 start1 end1 sample2 start2 end2 type
    Returns a DataFrame with calculated length column.
    """
    cols = ["sample1", "start1", "end1", "sample2", "start2", "end2", "type"]
    df = pd.read_csv(filepath, sep="\t", names=cols, comment="#", engine="python")

    # Ensure numeric columns
    for c in ["start1", "end1", "start2", "end2"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Compute SV length
    df["length"] = df.apply(
        lambda row: (row["end2"] - row["start2"]) if row["type"] == "I" else (row["end1"] - row["start1"]),
        axis=1
    )

    df["source_file"] = os.path.basename(filepath)
    return df


def main():
    args = arg_parser()

    # Collect all .bed, .bedpe, or .txt files
    bed_files = [
        os.path.join(args.bed_dir, f)
        for f in os.listdir(args.bed_dir)
        if f.endswith(".bed") or f.endswith(".bedpe") or f.endswith(".txt")
    ]

    if not bed_files:
        print(f"No BED files found in {args.bed_dir}")
        return

    all_dfs = []
    for fp in bed_files:
        print(f"Reading {fp} ...")
        df = read_bed_file(fp)
        all_dfs.append(df)

    combined = pd.concat(all_dfs, ignore_index=True)

    # Filter out invalid or zero-length entries
    combined = combined[combined["length"].notnull() & (combined["length"] > 0)]

    print(f"Total SVs: {len(combined)} ({combined['type'].value_counts().to_dict()})")

    # Split by SV type
    insertions = combined.loc[combined["type"] == "I", "length"]
    deletions = combined.loc[combined["type"] == "D", "length"]

    # Prepare data for violin plot
    data = [insertions, deletions]
    labels = ["Insertions (I)", "Deletions (D)"]

    # Create the violin plot
    plt.figure(figsize=(8, 6))
    parts = plt.violinplot(
        data,
        showmeans=True,
        showextrema=True,
        showmedians=True
    )

    # Customize appearance
    colors = ["#66c2a5", "#fc8d62"]
    for pc, color in zip(parts['bodies'], colors):
        pc.set_facecolor(color)
        pc.set_edgecolor('black')
        pc.set_alpha(0.8)

    plt.yscale("log")
    plt.xticks([1, 2], labels)
    plt.ylabel("SV Length (bp, log scale)")
    plt.title("SV Length Distributions (Insertions vs Deletions)")
    plt.grid(axis='y', linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig(args.output, dpi=300)


if __name__ == "__main__":
    main()
