import argparse
import os
import pandas as pd
import re

def extract_chr_from_dirname(dirname):
    """
    Extract 'chrX' from directory name like msa_chrX_sim_cases_YYYYMMDD
    """
    match = re.search(r'msa_(chr[\w\d]+)_sim_cases_\d+', dirname)
    return match.group(1) if match else None

def find_chr_dirs(base_path):
    """
    Find all subdirectories matching pattern msa_{chr}_sim_cases_{date}/
    """
    chr_dirs = []
    for entry in os.listdir(base_path):
        full_path = os.path.join(base_path, entry)
        if os.path.isdir(full_path) and re.match(r'^msa_chr[\w\d]+_sim_cases_\d+$', entry):
            chr_dirs.append(full_path)
    return chr_dirs

def find_case_dirs(chr_dir):
    """
    Find case_* subdirectories inside a given msa_{chr}_sim_cases_* directory
    """
    case_dirs = []
    for entry in os.listdir(chr_dir):
        full_path = os.path.join(chr_dir, entry)
        if os.path.isdir(full_path) and re.match(r'^case_\d+$', entry):
            case_dirs.append(full_path)
    return sorted(case_dirs, key=lambda x: int(re.search(r'case_(\d+)', x).group(1)))

def read_all_tree_comparisons(base_dir):
    """
    Process all msa_{chr}_sim_cases_* directories and read tree_comparison.tsv files
    """
    all_dfs = []

    chr_dirs = find_chr_dirs(base_dir)

    for chr_dir in chr_dirs:
        chr_name = extract_chr_from_dirname(os.path.basename(chr_dir))
        if not chr_name:
            print(f"Warning: Could not extract chr from {chr_dir}")
            continue

        case_dirs = find_case_dirs(chr_dir)
        for case_dir in case_dirs:
            file_path = os.path.join(case_dir, 'tree_comparison.tsv')
            if os.path.isfile(file_path):
                df = pd.read_csv(file_path, sep='\t', header=None)
                df.columns = ['height', 'size', 'match']
                df['case'] = os.path.basename(case_dir)
                df['chr'] = chr_name
                all_dfs.append(df)
            else:
                print(f"Warning: Missing file {file_path}")

    return all_dfs

def main():
    parser = argparse.ArgumentParser(description="Parse tree_comparison.tsv files from simulation directories.")
    parser.add_argument('path', help='Path to simulations/MSA_simulations/')
    args = parser.parse_args()

    all_dfs = read_all_tree_comparisons(args.path)

    if not all_dfs:
        print("No tree_comparison.tsv files found.")
        return

    combined_df = pd.concat(all_dfs, ignore_index=True)
    print(combined_df.head())

    # Optional: save to file
    # combined_df.to_csv("all_tree_comparisons.tsv", sep='\t', index=False)

if __name__ == '__main__':
    main()
