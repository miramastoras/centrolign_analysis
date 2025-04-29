from ete3 import Tree, TreeStyle
import numpy as np
import itertools
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
import argparse
import csv

def arg_parser():
    '''
    Parses command line arguments with argparse
    '''
    parser = argparse.ArgumentParser(
        prog='regression_tree_eval.py',
        description="""Linear regression to determine tree goodness of fit for pairwise distances. Uses Jordan's proposal of 0 1 encoding of pairs in subtree for each branch""")

    parser.add_argument("-t", "--tree",
                        required=True,
                        help="tree in newick format")
    parser.add_argument("-s", "--samples",
                        required=True,
                        help="list of samples to plot from tree")
    parser.add_argument("-d", "--pairwise_distances",
                        required=True,
                        help="csv of sample1,sample1,dist")
    parser.add_argument("-f", "--distance_format",
                        required=False,
                        default="col",
                        help="either")
    parser.add_argument("-o", "--output_dir",
                        required=True,
                        help="directory path to write output files to")

    return parser.parse_args()

args = arg_parser()

# Prepare test tree
tree = Tree(args.tree)
with open(args.samples, 'r') as file:
    samples = [line.strip() for line in file]
tree.prune(samples, preserve_branch_length=True)

# assign names to nodes based on level order traversal for debugging
for i, node in enumerate(tree.traverse("levelorder")):
    if not node.name:  # If the node has no name
        node.name = f"Node_{i}"  # Assign a new name, e.g., "Node_0", "Node_1", etc.
    #print(f"Node name: {node.name}")

# Initialize counters
internal_node_count = 0
leaf_count = 0

# Traverse the tree and count internal nodes and leaf nodes
for node in tree.traverse():
    if node.is_leaf():
        leaf_count += 1  # Count leaf nodes
    else:
        internal_node_count += 1  # Count internal nodes (non-leaves)

# Print the results
print(f"Number of internal nodes: {internal_node_count}")
print(f"Number of leaves: {leaf_count}")

### X Matrix #####

# initialize matrix of zeroes
num_pairwise_combinations = leaf_count * (leaf_count - 1) // 2
x_matrix = np.zeros((num_pairwise_combinations, internal_node_count))

# get list of all pairwise combinations of sample names
leaf_names = [leaf.name for leaf in tree.get_leaves()]
pairwise_combinations = list(itertools.combinations(leaf_names, 2))
sorted_pairwise_combinations = [tuple(sorted(combination)) for combination in pairwise_combinations]

# get dictionary containing the index of each pairwise combination in the list
index_dict = {samples: index for index, samples in enumerate(sorted_pairwise_combinations)}

# for each internal node, build out vector of all pairwise combinations
col_index = 0
for node in tree.traverse("levelorder"):
    vector= np.zeros(num_pairwise_combinations) # initialize vector to zero
    if not node.is_leaf():  # Check if the node is not a leaf (internal nodes only)
        #print(f"Internal Node: {node.name}")
        # get all descendants of internal node
        descendant_names = [descendant.name for descendant in node.get_descendants() if descendant.is_leaf()]
        # get pairwise combinations of all descendants, with each tuple sorted by alphabetical order
        pairwise_descendants=[tuple(sorted(combination)) for combination in list(itertools.combinations(descendant_names, 2))]

        for comb in pairwise_descendants:
            vector_index=index_dict[comb] # row in vector belonging to this combination
            vector[vector_index]=1 # set vector to 1 indicating this pair is descended from this internal node

        x_matrix[:, col_index]=vector
        col_index+=1

# spot check
# col_index = 0
# for node in tree.traverse("levelorder"):
#     if not node.is_leaf():
#         #print(f"Internal Node: {node.name}")
#         for n in range(0,num_pairwise_combinations):
#             #print( sorted_pairwise_combinations[n] , x_matrix[index_dict[sorted_pairwise_combinations[n]] , col_index])
#         col_index += 1

### Y Matrix #####

y_values={}
# Read in pairwise distances to be the y variable
if args.distance_format == "matrix":
    distances_df = pd.read_csv(args.pairwise_distances, sep='\t', index_col=0,header=None)
    distances_df=distances_df.drop(distances_df.columns[-1], axis=1)
    distances_df.columns=distances_df.index

    for pairwise_comb in sorted_pairwise_combinations:
        y_values[tuple(sorted([pairwise_comb[0],pairwise_comb[1]]))]=distances_df.loc[pairwise_comb[0],pairwise_comb[1]]

else:
    with open(args.pairwise_distances, "r") as csvfile:
        reader = csv.reader(csvfile)
        headers = next(reader)  # Skip the header row
        for row in reader:
            y_values[tuple(sorted([row[0],row[1]]))]= float(row[2])

# create empty y matrix for pairwise distances
y_matrix=np.zeros(num_pairwise_combinations)

for pairwise_comb in sorted_pairwise_combinations:
    # index to place pairwise value in y matrix
    index=index_dict[pairwise_comb]

    # get distance value from csv
    dist=y_values[tuple(sorted([pairwise_comb[0],pairwise_comb[1]]))]

    # set dist value in y matrix
    y_matrix[index]=dist


# Initialize the Linear Regression model
model = LinearRegression()

# Train the model on the training data
model.fit(x_matrix, y_matrix)

y_pred = model.predict(x_matrix)

rmse = np.sqrt( mean_squared_error(y_matrix, y_pred))

# Compute the residuals
residuals = y_matrix - y_pred

print(f"Root Mean Squared Error: {rmse}")

with open(args.output_dir +'residuals_pairwise.csv', 'w') as file:
    print("sample1,sample2,residual",file=file)
    for index in range(0,num_pairwise_combinations):
        print(sorted_pairwise_combinations[index][0],sorted_pairwise_combinations[index][1], residuals[index],sep=",",file=file)