import argparse
import pandas as pd

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Read a CSV file and process a specific column.")
    parser.add_argument("qc_csv_file", help="Path to the QC CSV file")
    parser.add_argument("r2_csv", help="path to release 2 assembly CSV file")
    args = parser.parse_args()

    qc = pd.read_csv(args.qc_csv_file)
    r2 = pd.read_csv(args.r2_csv)

    print(f"Rows: {qc.shape[0]}")
    print(f"Columns: {qc.shape[1]}")
    print(qc.columns)
    print(r2.columns)

    # Merge in s3 links
    qc_merged = qc.merge(r2[['assembly_name', 'assembly', 'assembly_fai']], on='assembly_name', how='left')

    print(f"Rows: {qc_merged.shape[0]}")
    print(f"Columns: {qc_merged.shape[1]}")


if __name__ == "__main__":
    main()
