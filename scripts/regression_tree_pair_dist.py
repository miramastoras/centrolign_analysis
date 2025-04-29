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
        description="""Linear regression to determine tree goodness of fit for pairwise distances. Uses Benedict's proposal with pairwise distance according to branch length""")

    parser.add_argument("-t", "--pair_tree_dist",
                        required=True,
                        help="pairwise distances from the tree according to tree_pair_dist centrolign script")
    parser.add_argument("-s", "--samples",
                        required=True,
                        help="list of samples to plot from tree")
    parser.add_argument("-p", "--all_pairs",
                        required=True,
                        help="csv of sample1,sample1,dist from centrolign all pairs distances")
    parser.add_argument("-o", "--output_dir",
                        required=True,
                        help="directory path to write output files to")

    return parser.parse_args()

args = arg_parser()

### X Matrix #####
x_values={}
# Read in pairwise distances to be the y variable
with open(args.pair_tree_dist, "r") as csvfile:
    reader = csv.reader(csvfile)
    headers = next(reader)  # Skip the header row
    for row in reader:
        x_values[tuple(sorted([row[0],row[1]]))]= float(row[2])

### Y Matrix #####
y_values={}
# Read in pairwise distances to be the y variable
with open(args.all_pairs, "r") as csvfile:
    reader = csv.reader(csvfile)
    headers = next(reader)  # Skip the header row
    for row in reader:
        y_values[tuple(sorted([row[0],row[1]]))]= float(row[2])

# use sample list to create list of pairwise combinations
with open(args.samples, 'r') as file:
    samples = [line.strip() for line in file]
pairwise_combinations = list(itertools.combinations(samples, 2))
sorted_pairwise_combinations = [tuple(sorted(combination)) for combination in pairwise_combinations]
num_pairwise_combinations = len(samples) * (len(samples) - 1) // 2

# get dictionary containing the index of each pairwise combination in the list
index_dict = {samples: index for index, samples in enumerate(sorted_pairwise_combinations)}

# create empty y matrix for pairwise distances
y_matrix=np.zeros(num_pairwise_combinations)
x_matrix=np.zeros(num_pairwise_combinations)

for pairwise_comb in sorted_pairwise_combinations:
    # index to place pairwise value in y matrix
    index=index_dict[pairwise_comb]

    # get distance value from csv
    dist=y_values[tuple(sorted([pairwise_comb[0],pairwise_comb[1]]))]
    # set dist value in y matrix
    y_matrix[index]=dist

    # get distance value from csv
    dist2 = x_values[tuple(sorted([pairwise_comb[0], pairwise_comb[1]]))]
    # set dist value in y matrix
    x_matrix[index] = dist2


# Reshape x to 2D for sklearn
x = x_matrix.reshape(-1, 1)

# Initialize the Linear Regression model
model = LinearRegression()

# Train the model on the training data
model.fit(x, y_matrix)

y_pred = model.predict(x)

rmse = np.sqrt( mean_squared_error(y_matrix, y_pred))

# Compute the residuals
residuals = y_matrix - y_pred

print(f"Root Mean Squared Error: {rmse}")

with open(args.output_dir +'residuals_pairwise.csv', 'w') as file:
    print("sample1,sample2,residual",file=file)
    for index in range(0,num_pairwise_combinations):
        print(sorted_pairwise_combinations[index][0],sorted_pairwise_combinations[index][1], residuals[index],sep=",",file=file)