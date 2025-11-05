#!/usr/bin/env python3

import argparse
import os
import pandas as pd

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
        "-o", "--output",
        required=False,
        default="merged_sv_beds.tsv",
        help="Output TSV file (default: merged_sv_beds.tsv)"
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
                all_beds.append(df)
            except Exception as e:
                print("Warning: Could not read {}: {}".format(bed_path, e))

    if all_beds:
        return pd.concat(all_beds, ignore_index=True)
    else:
        print("Warning: No .bed files found in {} for clade '{}'".format(bed_folder, clade))
        return pd.DataFrame(columns=["sample1", "start1", "end1", "sample2", "start2", "end2", "type", "diff", "clade", "source_file"])


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
        # Print head of DataFrame
        print("=== Head of merged DataFrame ===")
        print(merged_df.head(), "\n")

        # Print count of rows per clade
        print("=== Record count per clade ===")
        print(merged_df["clade"].value_counts())
    else:
        print("No data merged; check your input CSV and folder paths.")


if __name__ == "__main__":
    main()
