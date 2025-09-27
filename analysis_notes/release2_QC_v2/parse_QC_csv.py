import argparse
import pandas as pd

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Read a CSV file and process a specific column.")
    parser.add_argument("qc_csv_file", help="Path to the QC CSV file")
    parser.add_argument("r2_csv", help="path to release 2 assembly CSV file")
    parser.add_argument("new_csv_name", help="path to new csv file name")
    args = parser.parse_args()

    # read in Julian's QC csv and the R2 assembly CSV with s3 links
    qc = pd.read_csv(args.qc_csv_file)
    r2 = pd.read_csv(args.r2_csv)

    # Merge in s3 links
    qc_merged = qc.merge(r2[['assembly_id', 'assembly', 'assembly_fai']], on='assembly_id', how='left')

    # write only necessary columns to file
    columns_to_write = ['sample_id', 'haplotype', 'assembly_id', 'sequence_id', 'asat_start',
       'asat_end', 'assembly', 'assembly_fai','chrom_assignment']

    qc_merged.to_csv(args.new_csv_name, columns=columns_to_write, index=False)


if __name__ == "__main__":
    main()
