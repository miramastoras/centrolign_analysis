
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

plt.show()