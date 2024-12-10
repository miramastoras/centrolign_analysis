#!/usr/bin/env python3

'''
Purpose: Plot cenhap (flank) distance vs alignment distance
Author: Mira Mastoras, mmastora@ucsc.edu
python3 distance_compare.py \
        -c /Users/miramastoras/Desktop/distance_compare/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_hap.hD.txt \
        -a /Users/miramastoras/Desktop/distance_compare/pairwise_distance.csv \
        -o /Users/miramastoras/Desktop/distance_compare/
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
        prog='distance_compare.py',
        description="""Plot cenhap (flank) distance vs alignment distance""")

    parser.add_argument("-c", "--cenhap_dist",
                        required=True,
                        help="tab separated matrix of pairwise distance values from flanks")
    parser.add_argument("-a", "--alignment_dist",
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
    alignment_dists = {}
    with open(args.alignment_dist, "r") as csvfile:
        reader = csv.reader(csvfile)
        headers = next(reader)  # Skip the header row
        for row in reader:
            key = "_".join(sorted([row[0], row[1]]))
            value = float(row[2])  # Column 3 as the value
            alignment_dists[key] = value

    # Read in flank distance values
    df = pd.read_csv(args.cenhap_dist, sep='\t', index_col=0)

    # initialize plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 6),gridspec_kw = {'wspace':0, 'hspace':0,'width_ratios': [4, 1]})

    # initialize lists for getting max and min values
    cenhap_list=[]
    aln_list=[]

    num_keys=len(alignment_dists.keys())
    rainbow_colormap = cm.get_cmap('rainbow', num_keys+1)
    colors = [rainbow_colormap(i) for i in range(num_keys)]
    # Loop through pairs in alignment distance dict

    pair_num=0
    markers=['o','s','^','v','D']
    markers_dup=markers * num_keys

    axes[1].set_xlim(0,1)
    for key in alignment_dists.keys():
        pair1=str(key.split("_")[0])
        pair2=str(key.split("_")[1])

        if pair1 != pair2:
            aln=alignment_dists[key]
            aln_list.append(aln)
            cenhap=df.loc[pair1,pair2]
            cenhap_list.append(cenhap)

            axes[0].plot(cenhap,aln,markersize=3,color=colors[pair_num],marker=markers_dup[pair_num])
            axes[1].plot(.1,pair_num,markersize=5,color=colors[pair_num],marker=markers_dup[pair_num])
            axes[1].text(0.15,pair_num, pair1+","+pair2, fontsize=8, va="center")
            pair_num+=1

    axes[0].set_xlabel('Cenhap (flank) distance')
    axes[0].set_ylabel('Alignment based distance')
    axes[1].set_axis_off()
    plt.show()
    fig.savefig(args.output_dir + "alignment_vs_cenhap_dist.png", dpi=600)
if __name__ == '__main__':
    main()
