import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import toytree

def plot_tree_and_half_matrix_fast(newick, lower_data, upper_data, labels, figsize=(25, 20),
                                    cmap_lower='Blues', cmap_upper='Reds'):

    # Load the Newick tree
    tree = toytree.tree(newick)

    # Get the leaf order from the tree
    leaf_names = tree.get_tip_labels()

    # Map labels to indices
    label_to_idx = {label: idx for idx, label in enumerate(labels)}
    assert set(leaf_names) <= set(labels), "All tree leaves must be present in the labels."

    # Reorder matrices according to tree order
    indices = [label_to_idx[name] for name in leaf_names]
    lower_ordered = lower_data[np.ix_(indices, indices)]
    upper_ordered = upper_data[np.ix_(indices, indices)]

    # Create the plot
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(1, 5, width_ratios=[1.5, 0.2, 10, 0.2, 0.2], wspace=0.05)

    # Tree plot
    ax_tree = fig.add_subplot(gs[0])
    tree.draw(
        ax=ax_tree,
        use_edge_lengths=True,
        edge_widths=1,
        tip_labels_align=True,
        tip_labels=True,
        tip_labels_fontsize=6,
    )
    max_dist = tree.treenode.get_farthest_leaf_topology()[1]
    ax_tree.set_xlim(-0.1, max_dist + 0.5)
    ax_tree.set_ylim(-1, len(leaf_names))
    ax_tree.axis('off')

    # Colorbars placeholders
    ax_cbar_lower = fig.add_subplot(gs[1])
    ax_cbar_upper = fig.add_subplot(gs[4])

    # Matrix plot
    ax_matrix = fig.add_subplot(gs[2])

    mask_lower = np.triu(np.ones_like(lower_ordered, dtype=bool))
    mask_upper = np.tril(np.ones_like(upper_ordered, dtype=bool), -1)

    # Lower triangle
    norm_lower = plt.Normalize(vmin=np.nanmin(lower_ordered), vmax=np.nanmax(lower_ordered))
    sns.heatmap(lower_ordered, mask=mask_lower, cmap=cmap_lower, square=True,
                cbar_ax=ax_cbar_lower, cbar_kws={'orientation': 'vertical'},
                ax=ax_matrix, linewidths=0.1, linecolor='white', norm=norm_lower)

    # Upper triangle
    norm_upper = plt.Normalize(vmin=np.nanmin(upper_ordered), vmax=np.nanmax(upper_ordered))
    sns.heatmap(upper_ordered, mask=mask_upper, cmap=cmap_upper, square=True,
                cbar_ax=ax_cbar_upper, cbar_kws={'orientation': 'vertical'},
                ax=ax_matrix, linewidths=0.1, linecolor='white', norm=norm_upper)

    n = len(leaf_names)
    ax_matrix.set_xticks(np.arange(n) + 0.5)
    ax_matrix.set_yticks(np.arange(n) + 0.5)
    ax_matrix.set_xticklabels(leaf_names, rotation=90, fontsize=5)
    ax_matrix.set_yticklabels(leaf_names, rotation=0, fontsize=5)

    # Manually adjust layout
    fig.subplots_adjust(left=0.05, right=0.97, top=0.97, bottom=0.05, wspace=0.1)
    plt.show()

# ----------------------------
# Example usage
np.random.seed(42)
n_samples = 50  # Try 500 or 1000

labels = [f"Taxon_{i}" for i in range(n_samples)]

lower = np.random.rand(n_samples, n_samples)
upper = np.random.rand(n_samples, n_samples)
lower = (lower + lower.T) / 2
upper = (upper + upper.T) / 2

# Random Newick tree
newick = "(" + ",".join(labels) + ");"

plot_tree_and_half_matrix_fast(newick, lower, upper, labels)