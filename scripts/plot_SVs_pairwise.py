#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def arg_parser():
    '''
    Parses command line arguments with argparse
    '''
    parser = argparse.ArgumentParser(
        description="Combine pairwise SV BED files by clade into a single dataframe."
    )
    parser.add_argument(
        "-i", "--input_csv",
        required=True,
        help="CSV file with columns: clade,path_to_SV_beds"
    )
    parser.add_argument(
        "-p", "--plot_prefix",
        required=False,
        default="sv_length_dist",
        help="Output prefix for generated plot PNG files (default: sv_length_dist)"
    )
    return parser.parse_args()

def read_sv_bed_files(clade, bed_folder):
    """
    Read all .bed files in a folder and add a 'clade' column.
    Each bed file has columns:
        sample1, start, end, sample2, start, end, type, diff
    """
    all_beds = []

    # Ensure folder exists
    if not os.path.isdir(bed_folder):
        raise FileNotFoundError("Folder not found for clade '{}': {}".format(clade, bed_folder))

    for filename in os.listdir(bed_folder):
        if filename.endswith(".bed"):
            bed_path = os.path.join(bed_folder, filename)
            try:
                df = pd.read_csv(
                    bed_path,
                    sep="\t",
                    header=None,
                    names=["sample1", "start1", "end1", "sample2", "start2", "end2", "type", "diff"]
                )
                df["clade"] = clade
                df["source_file"] = filename  # optional, helps track origin
                df["length"] = np.where(
                    df["type"] == "I",
                    df["end2"] - df["start2"],  # insertion → use sample2 coords
                    df["end1"] - df["start1"]  # otherwise (deletion) → sample1 coords
                )
                all_beds.append(df)
            except Exception as e:
                print("Warning: Could not read {}: {}".format(bed_path, e))

    if all_beds:
        return pd.concat(all_beds, ignore_index=True)
    else:
        print("Warning: No .bed files found in {} for clade '{}'".format(bed_folder, clade))
        return pd.DataFrame(columns=["sample1", "start1", "end1", "sample2", "start2", "end2", "type", "diff", "clade", "source_file"])

def plot_length_histogram(df, category_label, cond, output_file, bin_size=1000, max_bin=50000):
    """
    Plot SV length histogram for a single category.
    Bars aligned to bin start; last bin = '>1Mb'.
    """
    lengths = df[cond]["length"].values
    if lengths.size == 0:
        print(f"⚠️ No data for {category_label}, skipping plot.")
        return

    # Cap lengths
    lengths_capped = np.where(lengths >= max_bin, max_bin, lengths)

    # Compute histogram
    bins = np.arange(0, max_bin + bin_size, bin_size)
    counts, _ = np.histogram(lengths_capped, bins=bins)
    overflow_count = (lengths >= max_bin).sum()
    counts[-1] += overflow_count

    # Print counts per bin
    print(f"\n=== {category_label} ===")
    for i, c in enumerate(counts[:-1]):
        print(f"{int(bins[i])} | {c}")
    print(f">=50kb | {overflow_count}")

    # Plot histogram
    plt.figure(figsize=(12, 5))
    plt.bar(bins[:-1], counts, width=bin_size, align='edge', color='skyblue', edgecolor='black', alpha=0.7)

    # Labels at start of each bar
    labels = [str(int(b)) for b in bins[:-1]] + [">50kb"]
    plt.xticks(list(bins[:-1]) + [max_bin], labels, rotation=90)
    plt.xlabel("Length (bp)")
    plt.ylabel("Count")
    plt.title(f"SV Length Histogram: {category_label}")
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print(f"✅ Saved histogram: {output_file}")

def plot_category_counts(df, output_file):
    """
    Plot total SV counts for six categories:
    x-axis = category
    y-axis = SV count
    """
    import matplotlib.pyplot as plt

    if df.empty:
        print("⚠️ DataFrame is empty, skipping plot.")
        return

    categories = [
        ("I_diff_eq_-1", (df["type"] == "I") & (df["diff"] == -1)),
        ("I_diff_lt_0.1", (df["type"] == "I") & (df["diff"] < 0.1)),
        ("I_diff_gt_0.1", (df["type"] == "I") & (df["diff"] > 0.1)),
        ("D_diff_eq_-1", (df["type"] == "D") & (df["diff"] == -1)),
        ("D_diff_lt_0.1", (df["type"] == "D") & (df["diff"] < 0.1)),
        ("D_diff_gt_0.1", (df["type"] == "D") & (df["diff"] > 0.1)),
    ]

    counts = []
    labels = []

    for label, cond in categories:
        count = df[cond].shape[0]
        counts.append(count)
        labels.append(label)
        print(f"{label}: {count} SVs")

    # Plot
    plt.figure(figsize=(8, 5))
    plt.bar(labels, counts, color="skyblue",align='edge', edgecolor="black", alpha=0.7)
    plt.xlabel("Category")
    plt.ylabel("SV count")
    plt.title("SV Counts by Category")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print(f"✅ Saved plot: {output_file}")


def main():
    # parse command line arguments
    args = arg_parser()

    # Read input CSV
    clade_info = pd.read_csv(args.input_csv)
    if not {"clade", "path_to_SV_beds"}.issubset(clade_info.columns):
        raise ValueError("Input CSV must contain columns: 'clade' and 'path_to_SV_beds'")

    all_data = []

    for _, row in clade_info.iterrows():
        clade = str(row["clade"])
        bed_folder = str(row["path_to_SV_beds"])
        df = read_sv_bed_files(clade, bed_folder)
        if not df.empty:
            all_data.append(df)

    if all_data:
        merged_df = pd.concat(all_data, ignore_index=True)

        print("=== Head of merged DataFrame ===")
        print(merged_df.head(), "\n")

        # Print count of rows per clade
        print("=== Record count per clade ===")
        print(merged_df["clade"].value_counts())

        plot_length_distributions(merged_df, args.plot_prefix)

        top5 = merged_df.sort_values(by="length", ascending=False).head(5)
        print("\n=== Top 5 largest SVs ===")
        print(top5[["sample1", "sample2", "clade", "length", "type", "diff"]])

        output_file = f"{args.plot_prefix}_category_counts.png"
        plot_category_counts(merged_df, output_file)

    else:
        print("No data merged; check your input CSV and folder paths.")


if __name__ == "__main__":
    main()
