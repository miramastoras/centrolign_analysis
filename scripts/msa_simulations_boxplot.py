import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
import seaborn as sns


# ---- STEP 1: Read in the file ----
file_path = '/Users/miramastoras/Desktop/github_repos/centrolign_analysis/analysis_notes/simulations/summary_tables/msa_all_chroms_sim_cases_20250402_aln_summary_tables.txt'

# ---- STEP 2: Define column names ----
column_names = [
    "case", "distance", "truth_matches", "truth_match_rate", "matches", "match_rate",
    "mismatches", "mismatch_rate", "recall", "precision", "chr"
]

df = pd.read_csv(file_path, delim_whitespace=True, names=column_names)

# ---- STEP 3: Calculate F1 score ----
df['F1'] = 2 * (df['precision'] * df['recall']) / (df['precision'] + df['recall'])
df['F1'] = df['F1'].fillna(0)
df = df[df['F1'] > 0]  # optional: exclude zeros

# ---- STEP 4: Sort chromosomes numerically ----
def chr_sort_key(chr_name):
    match = re.match(r'chr(\d+)', chr_name)
    return int(match.group(1)) if match else float('inf')

chromosomes = sorted(df['chr'].unique(), key=chr_sort_key)

# ---- STEP 5: Set up plotting ----
fig, ax = plt.subplots(figsize=(12, 6))

box_width = 0.8
x_positions = np.arange(len(chromosomes))

fig, ax = plt.subplots(figsize=(12, 6))

#---- STEP 6: Plot each chromosome ----
for i, chr_ in enumerate(chromosomes):
    chr_data = df[df['chr'] == chr_]['F1'].values
    bplots = ax.boxplot(
        [chr_data],  # Must be a sequence, even if 1 box
        positions=[x_positions[i]],
        widths=box_width * 0.9,
        patch_artist=True,
        manage_ticks=False,
        showfliers=False,
        boxprops={'linewidth': 1.0},
        whiskerprops={'linewidth': 1.0},
        capprops={'linewidth': 1.0},
        medianprops={'color': 'black', 'linewidth': 1.0}
    )

    # Set facecolor manually
    for patch in bplots['boxes']:
        patch.set_facecolor('#56B4E9')  # blue




# ---- STEP 8: Format plot ----
ax.set_xticks(x_positions)
ax.set_xticklabels(chromosomes)
ax.set_xlim(-1, len(chromosomes))
ax.set_ylabel('F1 Score')
ax.set_title('Centrolign MSA simulations')

# ---- STEP 9: Save to file ----
plt.tight_layout()

plt.savefig('/Users/miramastoras/Desktop/github_repos/centrolign_analysis/analysis_notes/simulations/figures/msa_simulations_boxplots.png', dpi=600, bbox_inches='tight')
plt.savefig('/Users/miramastoras/Desktop/github_repos/centrolign_analysis/analysis_notes/simulations/figures/msa_simulations_boxplots.svg', dpi=600, bbox_inches='tight')
