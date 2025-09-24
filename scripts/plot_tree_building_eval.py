import argparse
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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
    #
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

    for i, y_vals in enumerate(violin_data):
        if len(y_vals) > 0:
            mean = np.mean(y_vals)
            ax.plot(i + 1, mean, marker='o', color='black', markersize=5, zorder=3)

    # sns.swarmplot(
    #     data=summary_df,
    #     x='chr',
    #     y='percent_correct',
    #     size=3,  # dot size
    #     alpha=0.7
    # )

    # --- Step 5: Format axes ---
    ax.set_xticks(np.arange(1, len(chromosomes) + 1))
    ax.set_xticklabels(chromosomes)
    ax.set_ylabel('Percentage of correct bipartitions per tree')
    ax.set_xlabel('Chromosome')
    ax.set_title('Percentage of correct bipartitions per tree')
    ax.set_ylim(0, 100)

    plt.tight_layout()
    # Optional save
    plt.savefig('/private/groups/patenlab/mira/centrolign/simulations/tree_building/percent_correct_nodes_violin.png', dpi=300)
    #plt.savefig('percent_correct_nodes_swarmplot.svg')

    # Aggregate counts of correct and incorrect nodes by chromosome
    summary = combined_df.groupby(['chr', 'correct']).size().unstack(fill_value=0)

    # Ensure columns 0 and 1 exist (incorrect and correct)
    if 0 not in summary.columns:
        summary[0] = 0
    if 1 not in summary.columns:
        summary[1] = 0

    # Sort chromosomes (optional, customize if needed)
    def chr_sort_key(chr_name):
        import re
        match = re.match(r'chr(\d+)', chr_name)
        return int(match.group(1)) if match else float('inf')

    summary = summary.reindex(sorted(summary.index, key=chr_sort_key))

    # Plot
    fig, ax = plt.subplots(figsize=(12, 6))

    # Bottom is incorrect counts (match == 0)
    ax.bar(summary.index, summary[0], label='Incorrect nodes', color='tomato')

    # Stack correct counts (match == 1) on top
    ax.bar(summary.index, summary[1], bottom=summary[0], label='Correct nodes', color='seagreen')

    # Formatting
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('Number of bipartitions')
    ax.set_title('Total number correct vs incorrect bipartitions per chromosome')
    ax.legend()

    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()
    plt.savefig('/private/groups/patenlab/mira/centrolign/simulations/tree_building/percent_correct_nodes_barplot.png', dpi=300)

    # ---- Assumes your DataFrame is called `df` with columns: 'chr', 'match', 'height' ----

    # ---- Assumes df has columns: 'chr', 'match', 'height' ----

    # Sort chromosomes numerically
    def chr_sort_key(chr_name):
        match = re.match(r'chr(\d+)', chr_name)
        return int(match.group(1)) if match else float('inf')

    combined_df['chr'] = combined_df['chr'].astype(str)
    chromosomes = sorted(combined_df['chr'].unique(), key=chr_sort_key)
    combined_df['chr'] = pd.Categorical(combined_df['chr'], categories=chromosomes, ordered=True)

    # Create plot
    plt.figure(figsize=(14, 6))
    sns.swarmplot(
        data=combined_df,
        x='chr',
        y='height',
        hue='correct',
        palette={0: 'tomato', 1: 'seagreen'},
        dodge=True,  # separates hue values side by side
        size=2,  # dot size
        alpha=0.7
    )

    # sns.violinplot(
    #     data=combined_df,
    #     x='chr',
    #     y='height',
    #     hue='correct',
    #     palette={0: 'tomato', 1: 'seagreen'},
    #     dodge=True,  # separates hue values side by side
    #     size=3,  # dot size
    #     alpha=0.7,
    #     inner="box"
    # )


    # Labels and formatting
    plt.xlabel("Chromosome")
    plt.ylabel("Node Height")
    plt.title("Bipartition correctness by internal node height")
    plt.xticks(rotation=45)
    plt.ylim(bottom=0)
    plt.legend(title='Match', labels=['Incorrect', 'Correct'], loc='upper right')
    plt.tight_layout()

    plt.savefig('/private/groups/patenlab/mira/centrolign/simulations/tree_building/node_height_swarm.png', dpi=300)

    plt.figure(figsize=(14, 6))
    # sns.swarmplot(
    #     data=combined_df,
    #     x='chr',
    #     y='num_leaves',
    #     hue='correct',
    #     palette={0: 'tomato', 1: 'seagreen'},
    #     dodge=True,  # separates hue values side by side
    #     size=2,  # dot size
    #     alpha=0.7
    # )

    sns.violinplot(
        data=combined_df,
        x='chr',
        y='num_leaves',
        dodge=True,  # separates hue values side by side
        size=3,  # dot size
        alpha=0.7,
        inner="box"
    )

    # Labels and formatting
    plt.xlabel("Chromosome")
    plt.ylabel("Number of leaves in subtree")
    plt.title("Bipartition correctness by number of leaves in subtree")
    plt.xticks(rotation=45)
    plt.ylim(bottom=0)
    plt.legend(title='Match', labels=['Incorrect', 'Correct'], loc='upper right')
    plt.tight_layout()

    plt.savefig('/private/groups/patenlab/mira/centrolign/simulations/tree_building/num_leaves_swarm.png', dpi=300)


if __name__ == '__main__':
    main()
