import argparse
import os
import pandas as pd
import re

def find_case_dirs(base_path):
    """
    Find all subdirectories matching case_{num}
    """
    case_dirs = []
    for entry in os.listdir(base_path):
        full_path = os.path.join(base_path, entry)
        if os.path.isdir(full_path) and re.match(r'^case_\d+$', entry):
            case_dirs.append(full_path)
    return sorted(case_dirs, key=lambda x: int(re.search(r'case_(\d+)', x).group(1)))

def extract_chr_from_path(path):
    """
    Extract chr identifier (e.g., chr7) from directory like msa_chr7_sim_cases_20250402
    """
    basename = os.path.basename(os.path.normpath(path))
    match = re.search(r'(chr[\w\d]+)', basename)
    return match.group(1) if match else 'unknown_chr'

def read_tree_comparisons(base_path):
    """
    For each case directory, read tree_comparison.tsv and annotate with case and chr
    """
    chr_name = extract_chr_from_path(base_path)
    case_dirs = find_case_dirs(base_path)
    all_dfs = []

    for case_dir in case_dirs:
        file_path = os.path.join(case_dir, 'tree_comparison.tsv')
        if os.path.isfile(file_path):
            df = pd.read_csv(file_path, sep='\t', header=None)
            df.columns = ['height', 'size', 'match']
            df['case'] = os.path.basename(case_dir)
            df['chr'] = chr_name
            all_dfs.append(df)
        else:
            print(f'Warning: {file_path} not found.')

    return all_dfs

def main():
    parser = argparse.ArgumentParser(description='Parse tree_comparison.tsv files from case_* directories.')
    parser.add_argument('path', help='Path to msa_{chr}_sim_cases_*/ directory containing case_* subdirectories.')
    args = parser.parse_args()

    all_dfs = read_tree_comparisons(args.path)

    if not all_dfs:
        print("No valid tree_comparison.tsv files found.")
        return

    combined_df = pd.concat(all_dfs, ignore_index=True)
    print(combined_df.head())

    # Optional: Save to file
    # combined_df.to_csv("combined_tree_comparisons.tsv", sep='\t', index=False)

if __name__ == '__main__':
    main()
