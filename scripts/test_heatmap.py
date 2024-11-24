
import argparse
import matplotlib.pyplot as plt
from ete3 import Tree, TreeStyle
from Bio import Phylo
import numpy as np
from scipy.cluster.hierarchy import dendrogram
from matplotlib.patches import Rectangle
import matplotlib.patches as patches
import csv

fig, axes = plt.subplots(1, 1, figsize=(12, 6))
axes.set_ylim([-1,10])
axes.set_xlim([-1,10])

diagonal=1
numbers = list(range(7))

positions = numbers

for pos1 in range(len(positions)):
    for pos2 in range(pos1 + 1, len(positions)):
        bottom=[((diagonal/2)*(pos2-pos1-1)),((pos1+pos2)/2)-(diagonal/2)]
        top=[bottom[0],bottom[1]+diagonal]
        left=[(bottom[0]-(diagonal/2)),(bottom[1]+(diagonal/2))]
        right=[(bottom[0]+(diagonal/2)),(bottom[1]+(diagonal/2))]

        # draw diamond for current pair
        diamond_xy=[bottom,left,top,right]
        diamond = patches.Polygon(diamond_xy, fill=True, edgecolor='red')
        axes.add_patch(diamond)

#plt.show()

from itertools import combinations
import random

samples = [
    "NA06985.1",
    "NA06991.1",
    "NA06997.2",
    "NA18501.1",
    "NA19103.1",
    "NA18498.2",
    "NA19983.2",
    "HG00408.2",
    "HG00447.1",
    "HG00464.1",
    "NA20128.1"
]

# Generate all pairs with random values
pairs = []
for pair in combinations(samples, 2):
    random_value = random.random()  # generates number between 0 and 1
    pairs.append((pair[0], pair[1], random_value))

# Print all pairs
for pair in pairs:
    print(f"{pair[0]},{pair[1]},{pair[2]:.3f}")