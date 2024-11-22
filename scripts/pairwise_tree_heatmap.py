#!/usr/bin/env python3

'''
Purpose: Plot tree with heatmap displaying pairwise input values
Author: Mira Mastoras, mmastora@ucsc.edu
python3 pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/tree_heatmap/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.nwk.txt \
        -s /Users/miramastoras/Desktop/tree_heatmap/samples.txt
'''

import argparse
import matplotlib.pyplot as plt
from ete3 import Tree, TreeStyle
from Bio import Phylo
import numpy as np
from scipy.cluster.hierarchy import dendrogram
from matplotlib.patches import Rectangle

def arg_parser():
    '''
    Parses command line arguments with argparse
    '''
    parser = argparse.ArgumentParser(
        prog='pairwise_tree_heatmap.py',
        description="""Plot tree with heatmap displaying pairwise input values""")

    parser.add_argument("-t", "--tree",
                        required=True,
                        help="tree in newick format")
    parser.add_argument("-s", "--samples",
                        required=True,
                        help="list of samples to plot from tree")

    return parser.parse_args()

def tree_to_linkage_matrix(biopython_tree):
    # definition of linkage matrix:
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
    # code from stack overflow https://stackoverflow.com/questions/31033835/newick-tree-representation-to-scipy-cluster-hierarchy-linkage-matrix-format

    # create tree object via Biopython (Bio.Phylo)
    tree = biopython_tree

    # calculate tree height
    tree_height = max(tree.distance(c) for c in tree.find_clades(terminal=True))
    # map labels used in linkage matrix to sample names:
    id_map = {}
    for i, c in enumerate(tree.find_clades(terminal=True)):
        c.comment = (i, 1)
        id_map[i] = c.name
    print(id_map)
    #Phylo.draw(tree)
    # Internal nodes (ancestors) are collected with their distances from the root and sorted in descending order of distance.
    anc_lst = []
    for c in tree.find_clades(terminal=False):
        d = tree.distance(c)
        anc_lst.append((c, list(c), d))
    anc_lst.sort(key=lambda x: x[2], reverse=True)
    print(anc_lst)
    # running number of node
    nodes = len(list(tree.find_clades(terminal=True)))
    lnk_lst = []
    for anc, children, anc_d in anc_lst:
        n_children = len(children)
        assert n_children >= 2
        child1 = children[0]
        for child2 in children[1:]:
            id1, n_leaves1 = child1.comment
            id2, n_leaves2 = child2.comment
            total_leaves = n_leaves1 + n_leaves2
            anc.comment = (nodes, total_leaves)
            distance = tree_height - anc_d
            nodes += 1
            row = [id1, id2, distance, total_leaves]
            lnk_lst.append(row)
            child1 = anc
    return np.array(lnk_lst), id_map

def main():
    # parse command line arguments
    args = arg_parser()

    # read in samples we want to plot
    with open(args.samples, 'r') as file:
        samples = [line.strip() for line in file]

    # Prune the tree to retain only the target samples
    tree = Tree(args.tree, format=1)
    tree.prune(samples, preserve_branch_length=True)

    # write out pruned tree to file
    with open("/Users/miramastoras/Desktop/pruned_tree.nwk", "w") as output_file:
        output_file.write(tree.write())

    # read back in pruned tree
    tree = Phylo.read("/Users/miramastoras/Desktop/pruned_tree.nwk", "newick")

    # Plot two matplotlib grids side by side
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw = {'wspace':0, 'hspace':0})

    # convert newick to linkage matrix
    linkage_matrix, id_map = tree_to_linkage_matrix(tree)

    # plot dendogram on left axes
    dn1 = dendrogram(linkage_matrix, ax=axes[0], orientation='left',labels=[id_map.get(i, str(i)) for i in range(len(id_map))])

    axes[0].set_title('Pruned Guide Tree')
    axes[1].set_title('Pairwise Metric')

    # Remove x and y axis labels
    axes[1].set_xlabel('')  # Remove x axis label
    axes[1].set_ylabel('')  # Remove y axis label
    #axes[1].set_xticks([])  # Remove x ticks
    axes[1].set_yticks([])  # Remove y ticks

    leaf_label_y_map = {}
    for tick in axes[0].get_yticklabels():
        tick.set_verticalalignment('bottom')
        leaf_label_y_map[str(tick.get_text())] = int(tick.get_position()[1])
    # get mapping of y position of leaves and leaf labels


    print(leaf_label_y_map)

    # plot labels
    for sample in leaf_label_y_map.keys():
        plt.plot([0,0.1], [leaf_label_y_map[sample],leaf_label_y_map[sample]], linestyle=':')
        axes[1].text(0, leaf_label_y_map[sample], sample, fontsize=6, color='black',va='bottom')

    # set heatmap grid y axes to same scale as the tree grid, so labels match up.
    axes[1].set_ylim(float(axes[0].get_ylim()[0]), float(axes[0].get_ylim()[1]))
    axes[1].set_xlim(0,1)
    # Show the plot
    plt.show()
if __name__ == '__main__':
    main()
