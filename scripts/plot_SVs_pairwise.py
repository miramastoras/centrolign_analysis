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
    Generate plots showing SV length distributions for given conditions:
      1) type = "I", diff = -1
      2) type = "I", diff < 0.1
      3) type = "I", diff > 0.1
      4) type = "D", diff = -1
      5) type = "D", diff < 0.1
      6) type = "D", diff > 0.1
    """
    if df.empty:
        print("DataFrame is empty, skipping plots.")
        return

    # Define plotting conditions
    conditions = [
        ("I_diff_eq_-1", (df["type"] == "I") & (df["diff"] == -1)),
        ("I_diff_lt_0.1", (df["type"] == "I") & (df["diff"] < 0.1)),
        ("I_diff_gt_0.1", (df["type"] == "I") & (df["diff"] > 0.1)),
        ("D_diff_eq_-1", (df["type"] == "D") & (df["diff"] == -1)),
        ("D_diff_lt_0.1", (df["type"] == "D") & (df["diff"] < 0.1)),
        ("D_diff_gt_0.1", (df["type"] == "D") & (df["diff"] > 0.1)),
    ]

    # Ensure directory exists for the output prefix
    out_dir = os.path.dirname(output_prefix)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Plot using matplotlib only
    for label, cond in conditions:
        subset = df[cond]
        if subset.empty:
            print("No records for condition:", label)
            continue

        plt.figure(figsize=(8, 5))
        plt.hist(subset["length"], bins=50, color="skyblue", edgecolor="black", alpha=0.7)
        plt.title("SV Length Distribution: {}".format(label))
        plt.xlabel("Length (bp)")
        plt.ylabel("Count")
        plt.grid(axis="y", linestyle="--", alpha=0.7)
        plt.tight_layout()

        output_file = "{}_{}.png".format(output_prefix, label)
        plt.savefig(output_file)
        plt.close()
        print("Saved plot:", output_file)


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
    else:
        print("No data merged; check your input CSV and folder paths.")


if __name__ == "__main__":
    main()
