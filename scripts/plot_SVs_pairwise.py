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

def plot_length_distributions(df, output_prefix):
    """
    Generate SV length histograms with:
      - bin size = 100 bp
      - all SVs >= 1,000,000 bp in a single bin labeled '>1Mb'
    Also prints counts per bin.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import os

    if df.empty:
        print("⚠️ DataFrame is empty, skipping plots.")
        return

    # Ensure directory exists
    out_dir = os.path.dirname(output_prefix)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Define filter conditions
    conditions = [
        ("I_diff_eq_-1", (df["type"] == "I") & (df["diff"] == -1)),
        ("I_diff_lt_0.1", (df["type"] == "I") & (df["diff"] < 0.1)),
        ("I_diff_gt_0.1", (df["type"] == "I") & (df["diff"] > 0.1)),
        ("D_diff_eq_-1", (df["type"] == "D") & (df["diff"] == -1)),
        ("D_diff_lt_0.1", (df["type"] == "D") & (df["diff"] < 0.1)),
        ("D_diff_gt_0.1", (df["type"] == "D") & (df["diff"] > 0.1)),
    ]

    bin_size = 1000
    max_bin = 50000  # anything >= 1Mb goes in final bin

    for label, cond in conditions:
        subset = df[cond]
        if subset.empty:
            print("⚠️ No records for condition:", label)
            continue

        # Cap lengths at max_bin for plotting/counting
        lengths = subset["length"].copy()
        lengths_capped = np.where(lengths >= max_bin, max_bin, lengths)

        # Define bins: 0,100,200,...,1,000,000
        bins = np.arange(0, max_bin + bin_size, bin_size)

        counts, bin_edges = np.histogram(lengths_capped, bins=bins)

        # Count SVs >=1Mb separately and add to final bin
        overflow_count = (lengths >= max_bin).sum()
        counts[-1] += overflow_count

        # Print bin counts
        print(f"\n=== {label} ===")
        print("Bin range (bp) | Count")
        for i in range(len(counts) - 1):
            print(f"{int(bin_edges[i])}-{int(bin_edges[i+1]-1)} | {counts[i]}")
        print(f">=50,000 | {overflow_count}")

        # Plot histogram
        plt.figure(figsize=(10, 5))
        # x positions for bars
        x_pos = list(range(len(counts)))
        plt.bar(x_pos, counts, color="skyblue", edgecolor="black", alpha=0.7)
        # x-axis labels: bin start values and '>1Mb' for final bin
        x_labels = [str(int(bin_edges[i])) for i in range(len(counts) - 1)] + [">1Mb"]
        plt.xticks(x_pos, x_labels, rotation=90)
        plt.title(f"SV Length Distribution: {label}")
        plt.xlabel("Length (bp)")
        plt.ylabel("Count")
        plt.tight_layout()

        output_file = f"{output_prefix}_{label}.png"
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

    else:
        print("No data merged; check your input CSV and folder paths.")


if __name__ == "__main__":
    main()
