import argparse
import pandas as pd
import os

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Read a CSV file and process a specific column.")
    parser.add_argument("qc_csv_file", help="Path to the QC CSV file")
    parser.add_argument("output_path", help="Path to the output bed files")
    args = parser.parse_args()

    # read in Julian's QC csv
    df = pd.read_csv(args.qc_csv_file)

    # Group by assembly_id
    for assembly_id, group in df.groupby("assembly_id"):
        # Extract sample_id and haplotype (same for the whole group)
        sample_id = group["sample_id"].iloc[0]
        haplotype = group["haplotype"].iloc[0]

        # Prepare the output DataFrame
        bed_df = group[["sequence_id", "asat_start", "asat_end", "chrom_assignment"]]

        # Construct filename
        filename = f"{sample_id}.{haplotype}_asat_arrays.bed"

        # Write to BED file (tab-separated, no header/index)
        bed_df.to_csv(filename, sep="\t", header=False, index=False)


if __name__ == "__main__":
    main()
