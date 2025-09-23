import argparse
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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
                df.columns = ['height', 'num_leaves', 'correct']
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

    # --- Step 1: Aggregate ---
    summary_df = (
        combined_df.groupby(['chr', 'case'])
            .agg(total_nodes=('correct', 'count'), correct_matches=('correct', 'sum'))
            .reset_index()
    )
    summary_df['percent_correct'] = 100 * summary_df['correct_matches'] / summary_df['total_nodes']

    # --- Step 2: Sort chromosomes naturally ---
    def chr_sort_key(chr_name):
        match = re.match(r'chr(\d+)', chr_name)
        return int(match.group(1)) if match else float('inf')

    chromosomes = sorted(summary_df['chr'].unique(), key=chr_sort_key)

    # --- Step 3: Prepare data for violin plot ---
    violin_data = [summary_df[summary_df['chr'] == chr_]['percent_correct'].values for chr_ in chromosomes]

    # --- Step 4: Plot using Matplotlib ---
    fig, ax = plt.subplots(figsize=(12, 6))

    violin_parts = ax.violinplot(
        dataset=violin_data,
        showmeans=False,
        showmedians=True,  # Only show medians
        showextrema=False  # Hide min/max
    )

    # --- Optional: Custom styling for violins ---
    for vp in violin_parts['bodies']:
        vp.set_facecolor('#56B4E9')
        vp.set_edgecolor('black')
        vp.set_alpha(0.7)
        vp.set_linewidth(1)

    if 'cmedians' in violin_parts:
        violin_parts['cmedians'].set_color('black')
        violin_parts['cmedians'].set_linewidth(1.2)

    # --- Step 5: Format axes ---
    ax.set_xticks(np.arange(1, len(chromosomes) + 1))
    ax.set_xticklabels(chromosomes)
    ax.set_ylabel('Percent of Correct Nodes')
    ax.set_xlabel('Chromosome')
    ax.set_title('Distribution of Correct Internal Node Percentage per Case')
    ax.set_ylim(0, 100)

    plt.tight_layout()
    # Optional save
    plt.savefig('/private/groups/patenlab/mira/centrolign/simulations/tree_building/percent_correct_nodes_swarmplot.png', dpi=300)
    #plt.savefig('percent_correct_nodes_swarmplot.svg')

    #plt.show()
    # Optional: save to file
    # combined_df.to_csv("all_tree_comparisons.tsv", sep='\t', index=False)

if __name__ == '__main__':
    main()
