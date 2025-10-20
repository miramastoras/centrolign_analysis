import os
import re
import pandas as pd
import sys
import seaborn as sns
import matplotlib.pyplot as plt

def extract_chromosome(filename):
    match = re.search(r'(chr[0-9XYM]+)', filename)
    if match:
        return match.group(1)
    else:
        raise ValueError(f"Chromosome not found in filename: {filename}")

def read_pairwise_files(directory):
    all_rows = []
    for filename in os.listdir(directory):
        if filename.endswith('_pairwise_consistency.txt'):
            filepath = os.path.join(directory, filename)

            if os.path.getsize(filepath) == 0:
                print(f"Skipping empty file: {filename}")
                continue

            try:
                chrom = extract_chromosome(filename)
            except ValueError as e:
                print(f"Warning: {e}")
                continue

            try:
                df = pd.read_csv(filepath, sep="\t")
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

def plot_swarm(df, output_prefix):
    # Convert aligned_jaccard to numeric, coercing 'NA' or invalid entries to NaN
    df['aligned_jaccard'] = pd.to_numeric(df['aligned_jaccard'], errors='coerce')
    df = df.dropna(subset=['aligned_jaccard'])

    if df.empty:
        print("No valid aligned_jaccard values to plot.")
        return

    plt.figure(figsize=(12, 6))
    sns.set(style="whitegrid")
    ax = sns.swarmplot(data=df, x='chromosome', y='aligned_jaccard', size=.2)

    plt.xticks(rotation=45, ha='right')
    plt.title("Aligned Jaccard Similarity per Chromosome")
    plt.tight_layout()
    output_file = f"{output_prefix}_swarmplot.png"
    plt.savefig(output_file, dpi=300)
    print(f"Plot saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python combine_pairwise.py <directory_path> <output_prefix>")
        sys.exit(1)

    dir_path = sys.argv[1]
    output_prefix = sys.argv[2]

    combined_df = read_pairwise_files(dir_path)
    plot_swarm(combined_df, output_prefix)
