#!/usr/bin/env python3

'''
Purpose: Plot tree with heatmap displaying pairwise input values
Author: Mira Mastoras, mmastora@ucsc.edu
python3 pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/tree_heatmap/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.nwk.txt \
        -s /Users/miramastoras/Desktop/tree_heatmap/samples2.txt \
        -p /Users/miramastoras/Desktop/tree_heatmap/pairwise_combinations2.csv \
        -o 
'''

import argparse
import matplotlib.pyplot as plt
from ete3 import Tree, TreeStyle
from Bio import Phylo
import numpy as np
from scipy.cluster.hierarchy import dendrogram
from matplotlib.patches import Rectangle
import matplotlib.patches as patches
import csv

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
    parser.add_argument("-p", "--pairwise_values",
                        required=True,
                        help="csv file with pairwise sample combinations in col1-2 and heatmap value in col3")
    parser.add_argument("-m", "--metric_label",
                        required=False,
                        default="Alignment Distance",
                        help="label of pairwise metric to plot in heatmap")
    parser.add_argument("-o", "--output_dir",
                        required=True,
                        help="directory path to write output files to")

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

    #Phylo.draw(tree)
    # Internal nodes (ancestors) are collected with their distances from the root and sorted in descending order of distance.
    anc_lst = []
    for c in tree.find_clades(terminal=False):
        d = tree.distance(c)
        anc_lst.append((c, list(c), d))
    anc_lst.sort(key=lambda x: x[2], reverse=True)

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
    with open(args.output_dir+"pruned_tree.nwk", "w") as output_file:
        output_file.write(tree.write())

    # read back in pruned tree
    tree = Phylo.read(args.output_dir+"pruned_tree.nwk", "newick")

    pairwise_vals={}
    # Read in heatmap values
    with open(args.pairwise_values, "r") as csvfile:
        reader = csv.reader(csvfile)
        headers = next(reader)  # Skip the header row
        for row in reader:
            key = "_".join(sorted([row[0],row[1]]))
            value = float(row[2])  # Column 3 as the value
            pairwise_vals[key] = value

    min_pairwise = min(pairwise_vals.values())  # 10
    max_pairwise = max(pairwise_vals.values())

    # Plot two matplotlib grids side by side
    fig, axes = plt.subplots(1, 4, figsize=(12, 6), gridspec_kw = {'wspace':0, 'hspace':0,'width_ratios': [4, .5,2,.25]})

    # convert newick to linkage matrix
    linkage_matrix, id_map = tree_to_linkage_matrix(tree)

    # plot dendogram on left axes
    dn1 = dendrogram(linkage_matrix, ax=axes[0], orientation='left',labels=[id_map.get(i, str(i)) for i in range(len(id_map))],link_color_func=lambda x: 'black')

    # get mapping of y position of leaves and leaf labels
    leaf_label_y_map = {}
    for tick in axes[0].get_yticklabels():
        tick.set_verticalalignment('bottom') # required for label position to exactly match where leaf ends on dendogram
        leaf_label_y_map[str(tick.get_text())] = int(tick.get_position()[1])

    # the diamond diagonal is equal to the distance between adjacent leaves on the y axis
    diagonal = leaf_label_y_map[id_map[1]] - leaf_label_y_map[id_map[0]]
    offset = diagonal / 2

    # set heatmap grid y lim to same y lim as the tree grid, so labels match up. also the label grid in the middle
    # get y lim of grid 2 and 3 to match the dendogram, so everything lines up with leaves
    grid1_ylim=[float(axes[0].get_ylim()[0]), float(axes[0].get_ylim()[1])]

    axes[1].set_ylim(grid1_ylim)
    axes[2].set_ylim(grid1_ylim)

    # set x axis in grid 3 symmetrical so diamond is square, plus offset to show entire square for first column
    axes[2].set_xlim(float(axes[0].get_ylim()[0])-offset, float(axes[0].get_ylim()[1])/2)
    axes[1].set_xlim(0,1)
    axes[1].set_axis_off()

    # plot sample labels on axes 1
    label_ends=(float(axes[0].get_ylim()[1])) / 7
    for sample in leaf_label_y_map.keys():
        axes[1].plot([0,1], [leaf_label_y_map[sample],leaf_label_y_map[sample]], linestyle=':',color="black")
        axes[1].text(0, leaf_label_y_map[sample], sample, fontsize=7, color='black',va='bottom')

    # heatmap colors
    # Color scale
    seafoam = (159 / 255, 226 / 255, 191 / 255)
    deepblue = (0, 0, 139 / 255)
    R = np.linspace(seafoam[0], deepblue[0], 101)  # 101 steps going from 1 to 0
    G = np.linspace(seafoam[1], deepblue[1], 101)
    B = np.linspace(seafoam[2], deepblue[2], 101)

    # plot color scale in fourth grid
    for i in np.arange(0, 100, 1):
        rectangle = patches.Rectangle([0,( i / 100)], 1, 0.1,
                                         # move rectangle position with i , go in steps of 10 with 10 width
                                         facecolor=(R[i], G[i], B[i]),
                                         linewidth=0)
        axes[3].add_patch(rectangle)  # add rectangle to panel everytime we go through the loop

    axes[3].yaxis.tick_right()
    axes[3].yaxis.set_label_position("right")
    #axes[3].set_yticks([0,1])
    axes[3].set_xticklabels([])
    # normalize range of expression values between 0 and 100
    def normalizePoint(minVal, maxVal, x):
        return (int(((x - minVal) / (maxVal - minVal)) * 100))

    # loop through id_map position labels (0 to num samples - 1)
    positions = list(id_map.keys())

    # plot heatmap!
    for pos1 in range(len(positions)):
        for pos2 in range(pos1 + 1, len(positions)):
            dict_key="_".join(sorted([id_map[pos1],id_map[pos2]]))
            val=pairwise_vals[dict_key]
            norm_val = normalizePoint(min_pairwise, max_pairwise, float(val))

            y_pos1=leaf_label_y_map[id_map[pos1]]
            y_pos2=leaf_label_y_map[id_map[pos2]]
            bottom=[((diagonal/2)*(pos2-pos1-1)),((y_pos1+y_pos2)/2)-(diagonal/2)]
            top=[bottom[0],bottom[1]+diagonal]
            left=[(bottom[0]-(diagonal/2)),(bottom[1]+(diagonal/2))]
            right=[(bottom[0]+(diagonal/2)),(bottom[1]+(diagonal/2))]

            # draw diamond for current pair
            diamond_xy=[bottom,left,top,right]
            diamond = patches.Polygon(diamond_xy, fill=True, facecolor=(R[norm_val],G[norm_val],B[norm_val]),)
            axes[2].add_patch(diamond)

    # hide y axis labels for dendogram
    axes[0].set_yticklabels([])
    # Show the plot
    axes[2].set_axis_off()
    axes[0].yaxis.set_visible(False)  # Hide y-axis ticks and labels
    axes[0].spines['left'].set_visible(False)  # Hide the left spine
    axes[0].spines['right'].set_visible(False)  # Hide the right spine (optional)
    axes[0].spines['top'].set_visible(False)
    # Keep the x-axis
    axes[0].xaxis.set_visible(True)
    axes[0].set_title('Pruned Guide Tree')
    axes[2].set_title(args.metric_label)
    axes[0].set_xlabel('KYA')

    #plt.show()
    fig.savefig(args.output_dir+"pairwise_tree_heatmap.png", dpi=600)
if __name__ == '__main__':
    main()
