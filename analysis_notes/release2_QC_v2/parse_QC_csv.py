import argparse
import pandas as pd

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Read a CSV file and process a specific column.")
    parser.add_argument("qc_csv_file", help="Path to the QC CSV file")
    parser.add_argument("r2_csv", help="path to release 2 assembly CSV file")
    args = parser.parse_args()

    qc = pd.read_csv(args.qc_csv_file)
    r2 = pd.read_csv(args.qc_csv_file)

    print(qc.head())

    print(r2.head())

if __name__ == "__main__":
    main()
