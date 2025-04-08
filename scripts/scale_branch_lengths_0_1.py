from ete3 import Tree

tree = Tree("/Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj_upgma.nwk")

# A small epsilon to avoid zero branch lengths
epsilon = 1e-6

# Get all the branch lengths
branch_lengths = [n.dist for n in tree.traverse()]

# Shift the branch lengths by adding epsilon to avoid 0, and then rescale
min_length = min(branch_lengths)
max_length = max(branch_lengths)

# If all branch lengths are the same, scaling won't work. In that case, just return 1 for all.
if min_length == max_length:
    for node in tree.traverse():
        node.dist = 1
else:
    # Scale the branch lengths to be between epsilon and 1
    def scale_branch_length(length):
        # Add epsilon to avoid zero and rescale to the range (0, 1)
        return (length - min_length + epsilon) / (max_length - min_length + epsilon)

    # Apply the scaling function to each branch length
    for node in tree.traverse():
        node.dist = scale_branch_length(node.dist)

# Check the tree with scaled branch lengths
# tree.write(outfile="/Users/miramastoras/Desktop/Sasha_HPRC_trees_2_25_25/HPRC_chr12_P_Q_mira.2.25.25.upgma.rescaled_0_1.nwk", format=1)
tree.write(outfile="/Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj_upgma.rescaled_0_1.nwk", format=5)