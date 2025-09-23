import argparse
import os
import pandas as pd
import re

def find_case_dirs(base_path):
    """
    Return sorted list of subdirectories that match the pattern case_{num}
    """
    case_dirs = []
    for entry in os.listdir(base_path):
        full_path = os.path.join(base_path, entry)
        if os.path.isdir(full_path) and re.match(r'^case_\d+$', entry):
            case_dirs.append(full_path)
    return sorted(case_dirs, key=lambda x: int(re.search(r'case_(\d+)', x).group(1)))

def extract_chr_name(path):
    """
    Extract chromosome name like 'chr7' from the basename of the input path
    """
    basename = os.path.basename(os.path.normpath(path))
    match = re.search(r'(chr\d+)', basename)
    return match.group(1) if match else 'unknown_chr'

def read_tree_comparisons(base_path):
    """
    Find and read all tree_comparison.tsv files in case_* subdirectories
    """
    case_dirs = find_case_dirs(base_path)
    all_dfs = []

    chr_name = extract_chr_name(base_path)

    for case_dir in case_dirs:
        file_path = os.path.join(case_dir, 'tree_comparison.tsv')
        if os.path.isfile(file_path):
            df = pd.read_csv(file_path, sep='\t', header=None)
            df.columns = ['height', 'size', 'match']  # Adjust if needed
            df['case'] = os.path.basename(case_dir)
            df['chr'] = chr_name
            all_dfs.append(df)
        else:
            print(f'Warning: {file_path} not found.')

    return all_dfs

def main():
    parser = argparse.ArgumentParser(description='Read tree_comparison.tsv files from case_* subdirectories.')
    parser.add_argument('path', help='Base path containing case_* subdirectories')
    args = parser.parse_args()

    all_dfs = read_tree_comparisons(args.path)

    print(f"Found {len(all_dfs)} cases with valid tree_comparison.tsv files.")

    if all_dfs:
        combined_df = pd.concat(all_dfs, ignore_index=True)
        print(combined_df.head())
    else:
        print("No data found.")

if __name__ == '__main__':
    main()
