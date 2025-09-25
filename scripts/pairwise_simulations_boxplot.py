import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
import matplotlib.lines as mlines

# ---- STEP 1: Read in the file ----
# Update the file path if needed
file_path = '/Users/miramastoras/Desktop/github_repos/centrolign_analysis/analysis_notes/simulations/summary_tables/pair_all_chroms_sim_cases_20250421_aln_summary_tables.txt'

# ---- STEP 1: Define column names ----
column_names = [
    "case", "aligner", "distance", "truth_matches", "truth_match_rate",
    "matches", "match_rate", "mismatches", "mismatch_rate",
    "recall", "precision", "chr"
]

df = pd.read_csv(file_path, delim_whitespace=True, names=column_names)

# ---- STEP 2: Calculate F1 score ----
# F1 = 2 * (precision * recall) / (precision + recall)
# Make sure we handle division by zero properly
df['F1'] = 2 * (df['precision'] * df['recall']) / (df['precision'] + df['recall'])
df['F1'] = df['F1'].fillna(0)

# ----  Filter out F1 == 0 ----
df = df[df['F1'] > 0]

# ---- STEP 5: Get chr order and aligner groups ----
def chr_sort_key(chr_name):
    match = re.match(r'chr(\d+)', chr_name)
    return int(match.group(1)) if match else float('inf')

chromosomes = sorted(df['chr'].unique(), key=chr_sort_key)

aligners = sorted(df['aligner'].unique())
aligner_to_index = {aligner: i for i, aligner in enumerate(aligners)}

#  Compute mean F1 per (chr, aligner)
group_means = df.groupby(['chr', 'aligner'])['F1'].mean().reset_index()

# Set up boxplots
fig, ax = plt.subplots(figsize=(12, 6))

# Width of each boxplot group and each box
group_width = 0.8
box_width = group_width / len(aligners)

# Store box center positions to place stars
box_centers = {}

# For x-axis, positions of groups
x_group_positions = np.arange(len(chromosomes))

# Define colors for each aligner explicitly
aligner_colors = {
    'centrolign': '#56B4E9',  # blue
    'rama': '#E69F00',  # orange
    'unialigner': '#009E73',  # green
    # Add more if needed
}
# Even spacing around the center for each aligner
aligner_offsets = np.linspace(
    -group_width / 2 + box_width / 2,
    group_width / 2 - box_width / 2,
    len(aligners)
)

for i, chr_ in enumerate(chromosomes):
    chr_data = df[df['chr'] == chr_]

    data_to_plot = [chr_data[chr_data['aligner'] == aln]['F1'].values for aln in aligners]

    #positions = x_group_positions[i] - group_width/2 + box_width/2 + box_width * np.arange(len(aligners))
    positions = x_group_positions[i] + aligner_offsets

    bplots = ax.boxplot(
        data_to_plot,
        positions=positions,
        widths=box_width * 0.9,
        patch_artist=True,
        manage_ticks=False,
        showfliers=False,
        medianprops={'color': 'black'},
        boxprops={'linewidth': 1.0}  # <- Thinner box edge
        #whiskerprops={'linewidth': 0.8},  # <- Thinner whiskers
        #capprops={'linewidth': 0.8},  # <- Thinner caps
    )

    for patch, aln in zip(bplots['boxes'], aligners):
        patch.set_facecolor(aligner_colors.get(aln, 'gray'))  # fallback color
    # Save centers for stars
    for aln_pos, aln in zip(positions, aligners):
        box_centers[(chr_, aln)] = aln_pos

# Add stars to highlight one with highest mean
for chr_ in chromosomes:
    means_chr = group_means[group_means['chr'] == chr_]
    if means_chr.empty:
        continue

    best_row = means_chr.loc[means_chr['F1'].idxmax()]

    best_aligner = best_row['aligner']

    x_star = box_centers[(chr_, best_aligner)]
    max_f1 = df[(df['chr'] == chr_) & (df['aligner'] == best_aligner)]['F1'].max()
    y_star = max_f1 + 0.01  # smaller offset, closer to box top

    ax.text(
        x_star, y_star, '★',
        fontsize=10,  # smaller star size
        ha='center',
        va='bottom',
        color='black'
    )


# Formatting
ax.set_xticks(x_group_positions)
ax.set_xticklabels(chromosomes)
ax.set_xlim(-1, len(chromosomes))
ax.set_ylabel('F1 Score')
ax.set_title('Performance of simulated pairwise alignments')

# ---- Legend ----
legend_handles = [
    plt.Line2D([0], [0], color=aligner_colors.get(aln, 'gray'), lw=8)
    for aln in aligners
]

# Add star symbol to legend
star_handle = mlines.Line2D([], [], color='black', marker='*', linestyle='None', markersize=10, label='* = highest mean F1')

legend_handles.append(star_handle)

# Show legend at bottom
ax.legend(
    handles=legend_handles,
    labels=aligners + ['= highest mean F1'],
    loc='lower center',
    bbox_to_anchor=(0.5, -0.25),
    ncol=len(aligners) + 1,
    frameon=False
)

plt.tight_layout()
plt.savefig('/Users/miramastoras/Desktop/github_repos/centrolign_analysis/analysis_notes/simulations/figures/pairwise_simulations_boxplots.png', dpi=600, bbox_inches='tight')
plt.savefig('/Users/miramastoras/Desktop/github_repos/centrolign_analysis/analysis_notes/simulations/figures/pairwise_simulations_boxplots.svg', dpi=600, bbox_inches='tight')


# ===================== VIOLIN PLOT VERSION ==========================
fig_violin, ax_violin = plt.subplots(figsize=(12, 6))

# Recalculate centers for violin positions
violin_centers = {}

for i, chr_ in enumerate(chromosomes):
    chr_data = df[df['chr'] == chr_]

    data_to_plot = [chr_data[chr_data['aligner'] == aln]['F1'].values for aln in aligners]
    positions = x_group_positions[i] + aligner_offsets

    for j, (y_vals, aln) in enumerate(zip(data_to_plot, aligners)):
        parts = ax_violin.violinplot(
            dataset=y_vals,
            positions=[positions[j]],
            widths=box_width * 0.9,
            showmeans=False,
            showmedians=True,
            showextrema=False
        )
        # Set color for violin
        for pc in parts['bodies']:
            pc.set_facecolor(aligner_colors.get(aln, 'gray'))
            pc.set_alpha(0.8)
            pc.set_edgecolor('black')
            pc.set_linewidth(0.5)

        # Add manual median line
        median_val = np.median(y_vals)
        ax_violin.plot([positions[j] - box_width * 0.3, positions[j] + box_width * 0.3],
                       [median_val, median_val],
                       color='black', linewidth=1)

        # Save center for star
        violin_centers[(chr_, aln)] = positions[j]

# Add stars to highlight highest mean F1 per chromosome
for chr_ in chromosomes:
    means_chr = group_means[group_means['chr'] == chr_]
    if means_chr.empty:
        continue

    best_row = means_chr.loc[means_chr['F1'].idxmax()]
    best_aligner = best_row['aligner']

    x_star = violin_centers[(chr_, best_aligner)]
    max_f1 = df[(df['chr'] == chr_) & (df['aligner'] == best_aligner)]['F1'].max()
    y_star = max_f1 + 0.01

    ax_violin.text(
        x_star, y_star, '★',
        fontsize=10,
        ha='center',
        va='bottom',
        color='black'
    )

# Formatting
ax_violin.set_xticks(x_group_positions)
ax_violin.set_xticklabels(chromosomes)
ax_violin.set_xlim(-1, len(chromosomes))
ax_violin.set_ylabel('F1 Score')
ax_violin.set_title('Performance of simulated pairwise alignments')

# Legend (same as before)
legend_handles_violin = [
    plt.Line2D([0], [0], color=aligner_colors.get(aln, 'gray'), lw=8)
    for aln in aligners
]
star_handle_violin = mlines.Line2D([], [], color='black', marker='*', linestyle='None', markersize=10, label='* = highest mean F1')
legend_handles_violin.append(star_handle_violin)

ax_violin.legend(
    handles=legend_handles_violin,
    labels=aligners + ['= highest mean F1'],
    loc='lower center',
    bbox_to_anchor=(0.5, -0.25),
    ncol=len(aligners) + 1,
    frameon=False
)

plt.tight_layout()

# Save violin plot
plt.savefig('/Users/miramastoras/Desktop/github_repos/centrolign_analysis/analysis_notes/simulations/figures/pairwise_simulations_violins.png', dpi=600, bbox_inches='tight')
plt.savefig('/Users/miramastoras/Desktop/github_repos/centrolign_analysis/analysis_notes/simulations/figures/pairwise_simulations_violins.svg', dpi=600, bbox_inches='tight')
