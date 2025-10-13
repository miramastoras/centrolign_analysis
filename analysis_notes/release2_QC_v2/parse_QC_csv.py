import argparse
import pandas as pd
import os

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Reads concatenated asat csv files and outputs per sample beds for all the arrays")
    parser.add_argument("qc_table", help="Path to the QC table file. Either csv or tsv accepted")
    parser.add_argument("output_path", help="Path to the output bed files")
    args = parser.parse_args()

    # read in Julian's QC csv
    if args.qc_table.endswith(".tsv"):
        print("tsv")
        df = pd.read_csv(args.qc_table, sep="\t")

    elif args.qc_table.endswith(".csv"):
        df = pd.read_csv(args.qc_table)

    else:
        print("Error: invalid input file format. Must be either tsv or csv")
        exit()

    # Group by assembly_id
    for assembly_id, group in df.groupby("assembly_id"):
        # Extract sample_id and haplotype (same for the whole group)
        sample_id = group["sample_id"].iloc[0]
        haplotype = group["haplotype"].iloc[0]

        # Prepare the output DataFrame
        bed_df = group[["sequence_id", "asat_start", "asat_end", "chrom_assignment"]]

        # Construct filename
        filename = os.path.join(args.output_path, f"{sample_id}.{haplotype}_asat_arrays.bed")

        # Check if the file already exists
        file_exists = os.path.exists(filename)

        # Append if file exists, otherwise create new.
        # this is because chr 3 and chr4 are formatted differently
        bed_df.to_csv(
            filename,
            sep="\t",
            header=False,
            index=False,
            mode='a' if file_exists else 'w'
        )

if __name__ == "__main__":
    main()
