import os
import re
import pandas as pd
import sys

def extract_chromosome(filename):
    """
    Extract chromosome identifier from filename, e.g., 'chr12' from 'prefix_chr12_pairwise_consistency.txt'
    or 'prefix_chr12_group1_pairwise_consistency.txt'
    """
    match = re.search(r'(chr[0-9XYM]+)', filename)
    if match:
        return match.group(1)
    else:
        raise ValueError(f"Chromosome not found in filename: {filename}")

def read_pairwise_files(directory):
    """
    Find all *_pairwise_consistency.txt files in the directory, extract chromosome info,
    and load them into a single DataFrame. Skips empty files.
    """
    all_rows = []

    for filename in os.listdir(directory):
        if filename.endswith('_pairwise_consistency.txt'):
            filepath = os.path.join(directory, filename)

            # Skip if file is empty (0 bytes)
            if os.path.getsize(filepath) == 0:
                print(f"Skipping empty file: {filename}")
                continue

            try:
                chrom = extract_chromosome(filename)
            except ValueError as e:
                print(f"Warning: {e}")
                continue

            try:
                # Read the file
                df = pd.read_csv(filepath, sep="\t")

                # Skip if file has no rows
                if df.empty:
                    print(f"Skipping empty DataFrame: {filename}")
                    continue

                df['chromosome'] = chrom
                all_rows.append(df)

            except Exception as e:
                print(f"Error reading {filename}: {e}")
                continue

    if not all_rows:
        print("No valid non-empty files found.")
        return pd.DataFrame()

    return pd.concat(all_rows, ignore_index=True)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python plot_pairwise_consistency_per_chrom.py <directory_path>")
        sys.exit(1)

    dir_path = sys.argv[1]

    combined_df = read_pairwise_files(dir_path)

    # Print or save the result
    print(combined_df.tail())  # or save with combined_df.to_csv("output.csv", index=False)
