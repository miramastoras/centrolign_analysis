'''
Purpose: Compare centrolign pairwise distances from two experiments
Author: Mira Mastoras, mmastora@ucsc.edu
'''

import argparse
import matplotlib.pyplot as plt
import csv
import pandas as pd
import numpy as np
from matplotlib import cm

def arg_parser():
    '''
    Parses command line arguments with argparse
    '''
    parser = argparse.ArgumentParser(
        prog='compare_centrolign_pairwise_distances.py',
        description="""Compare centrolign pairwise distances from two experiments""")

    parser.add_argument("-a", "--alignment_dist1",
                        required=True,
                        help="csv file with pairwise distance values from centrolign in this format: Sample1,Sample2,Dist")
    parser.add_argument("-b", "--alignment_dist2",
                        required=True,
                        help="csv file with pairwise distance values from centrolign in this format: Sample1,Sample2,Dist")
    parser.add_argument("-o", "--output_dir",
                        required=True,
                        help="directory path to write output files to")

    return parser.parse_args()

def main():
    # parse command line arguments
    args = arg_parser()

    # Read in centrolign distance values
    dists_1 = {}
    with open(args.alignment_dist1, "r") as csvfile:
        reader = csv.reader(csvfile)
        headers = next(reader)  # Skip the header row
        for row in reader:
            key = "_".join(sorted([row[0], row[1]]))
            value = float(row[2])  # Column 3 as the value
            dists_1[key] = value

    dists_2 = {}
    with open(args.alignment_dist2, "r") as csvfile:
        reader = csv.reader(csvfile)
        headers = next(reader)  # Skip the header row
        for row in reader:
            key = "_".join(sorted([row[0], row[1]]))
            value = float(row[2])  # Column 3 as the value
            dists_2[key] = value

    # initialize plot
    fig, axes = plt.subplots(1, 1, figsize=(6, 6))

    # initialize lists for getting max and min values
    aln_list1=[]
    aln_list2=[]

    num_keys=len(dists_1.keys())
    rainbow_colormap = cm.get_cmap('rainbow', num_keys+1)
    colors = [rainbow_colormap(i) for i in range(num_keys)]
    # Loop through pairs in alignment distance dict

    pair_num=0
    markers=['o','s','^','v','D']
    markers_dup=markers * num_keys

    for key in dists_1.keys():
        pair1=str(key.split("_")[0])
        pair2=str(key.split("_")[1])

        if pair1 != pair2:
            aln1=dists_1[key]
            aln_list1.append(aln1)
            aln2=dists_2[key]
            aln_list2.append(aln2)

            axes.plot(aln1,aln2,markersize=1,color="black",marker="o")
            pair_num+=1

    plt.plot([0, 1], [0, 1],'--',color="red")
    axes.set_xlim([0,1.01])
    axes.set_ylim([0,1.01])
    axes.set_xlabel('Original Guide Tree')
    axes.set_ylabel('Tree infer method 1 (freezing deep nodes)')
    axes.set_title("Centrolign alignment distances")
    plt.show()
    fig.savefig(args.output_dir + "compare_centrolign_pairwise_dists.png", dpi=600)
if __name__ == '__main__':
    main()
