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
    and load them into a single DataFrame.
    """
    all_rows = []

    for filename in os.listdir(directory):
        if filename.endswith('_pairwise_consistency.txt'):
            filepath = os.path.join(directory, filename)
            try:
                chrom = extract_chromosome(filename)
            except ValueError as e:
                print(f"Warning: {e}")
                continue

            # Read the file into a DataFrame
            df = pd.read_csv(filepath, sep="\t")

            # Add chromosome column
            df['chromosome'] = chrom

            all_rows.append(df)

    if not all_rows:
        print("No valid files found.")
        return pd.DataFrame()

    # Concatenate all dataframes
    combined_df = pd.concat(all_rows, ignore_index=True)
    return combined_df

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python plot_pairwise_consistency_per_chrom.py <directory_path>")
        sys.exit(1)

    dir_path = sys.argv[1]

    combined_df = read_pairwise_files(dir_path)

    # Print or save the result
    print(combined_df.head())  # or save with combined_df.to_csv("output.csv", index=False)
