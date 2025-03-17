from ete3 import Tree,TreeStyle
import numpy as np
import itertools
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
import csv

# Prepare test tree
tree = Tree("/Users/miramastoras/Desktop/initial_experiments_local_computer/tree_heatmap/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.nwk.txt", format=5)
# Convert the DataFrame into a list of lists (each row as a list)
with open('/Users/miramastoras/Desktop/initial_experiments_local_computer/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.txt', 'r') as file:
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

# Read in pairwise distances to be the y variable

# sasha's format
#distances_df = pd.read_csv("/Users/miramastoras/Desktop/regression_trees/HPRC_chr12_P_Q_mira.2.25.25.m.tsv", sep='\t', index_col=0,header=None)
#distances_df=distances_df.drop(distances_df.columns[-1], axis=1)
#distances_df.columns=distances_df.index

pairwise_vals={}
# Read in heatmap values
with open("/Users/miramastoras/Desktop/initial_experiments_local_computer/centrolign_all_pairs/pairwise_distance.csv", "r") as csvfile:
    reader = csv.reader(csvfile)
    headers = next(reader)  # Skip the header row
    for row in reader:
        key = tuple(sorted([row[0],row[1]]))
        value = float(row[2])  # Column 3 as the value
        pairwise_vals[key] = value

# create empty y matrix for pairwise distances
y_matrix=np.zeros(num_pairwise_combinations)

for pairwise_comb in sorted_pairwise_combinations:
    # index to place pairwise value in y matrix
    index=index_dict[pairwise_comb]

    # get distance value from csv
    #dist=(distances_df.loc[pairwise_comb[0],pairwise_comb[1]]) # sasha's matrix format
    dist=pairwise_vals[pairwise_comb]
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

with open('/Users/miramastoras/Desktop/regression_trees/residuals_pairwise_chr12_old_tree_centrolign.csv', 'w') as file:
    print("sample1,sample2,residual",file=file)
    for index in range(0,num_pairwise_combinations):
        print(sorted_pairwise_combinations[index][0],sorted_pairwise_combinations[index][1], residuals[index],sep=",",file=file)