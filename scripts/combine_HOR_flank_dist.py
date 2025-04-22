'''
Purpose: Compute new pairwise distance metric, combining HOR distances derived from centrolign,
and flank SNP derived distances in a weighted sum
Author: Mira Mastoras, mmastora@ucsc.edu
'''

import argparse
import csv
import pandas as pd
import skbio

def arg_parser():
    '''
    Parses command line arguments with argparse
    '''
    parser = argparse.ArgumentParser(
        prog='combine_HOR_flank_dist.py',
        description="""Compute weighted sum of HOR and flank distances, infer neighbor joining tree from them""")

    parser.add_argument("-c", "--centrolign_HOR_dists",
                        required=True,
                        help="csv file with pairwise distance values from centrolign in this format: Sample1,Sample2,Dist")
    parser.add_argument("-f", "--flank_dists",
                        required=True,
                        help="comma separated matrix of pairwise distance values from flanks")
    parser.add_argument("-s", "--samples",
                        required=False,
                        help="txt file with samples to use.Default uses samples from centrolign HOR list.")
    parser.add_argument("-o", "--output_pre",
                        required=True,
                        help="directory path and file prefix to write output file to")

    return parser.parse_args()


def min_max_scale_df(df):
    """Scales all values in the DataFrame to be between 0 and 1, using global min and max."""
    df_scaled = df.copy().astype(float)  # Ensure numeric operations

    df_min = df.min().min()  # Global min across all values
    df_max = df.max().max()  # Global max across all values

    # Apply min-max scaling
    df_scaled = (df - df_min) / (df_max - df_min)

    return df_scaled

def main():

    # parse command line arguments
    args = arg_parser()

    # Read in centrolign distance values
    alignment_dists = {}

    # collect list of samples
    sample_list =[]

    with open(args.centrolign_HOR_dists, "r") as csvfile:
        reader = csv.reader(csvfile)
        headers = next(reader)  # Skip the header row
        for row in reader:
            sample_list.append(row[0])

            sample_list.append(row[1])
            key = "_".join(sorted([row[0], row[1]]))
            value = float(row[2])  # Column 3 as the value
            alignment_dists[key] = value
    print(type(sample_list[1]))
    if args.samples is not None:
        with open(args.samples, 'r') as file:
            samples = [line.strip() for line in file]
        sample_list=samples

    samps = sorted(set(sample_list))

    # Read in flank distance values
    flank_df = pd.read_csv(args.flank_dists, sep=',', index_col=0)

    # scale flank distance value to be between 0 and 1
    scaled_flank_df=min_max_scale_df(flank_df)

    # compute weighted sum, write to file, store in matrix for skbio
    mat = {}
    with open(args.output_pre +"_HOR_flank_dist_weighted.txt", "a") as file:
        for key in alignment_dists.keys():

            h=alignment_dists[key] # hor distance for each pairwise combination

            sample1=key.split("_")[0]
            sample2=key.split("_")[1]

            if sample1 in samps and sample2 in samps:
                f=scaled_flank_df.loc[sample1, sample2] # flank distance

                d = (1 - (1 - h)**2 + f**2) / 2

                # write to new file
                print(sample1, sample2, d, sep=",", file=file)

                # add to dictionary
                mat[(sample1, sample2)] = d
                mat[(sample2, sample1)] = d
            else:
                continue
    # reorganize as an array
    print(len(mat.keys()))
    D = []
    for samp1 in samps:
        D.append([])
        for samp2 in samps:
            if samp1 == samp2:
                D[-1].append(0.0)
            else:
                D[-1].append(mat[(samp1, samp2)])

    # make skbio type
    dist_mat = skbio.DistanceMatrix(D, samps)
    # print(dist_mat.to_data_frame(), file = sys.stderr)

    tree = skbio.tree.nj(dist_mat)
    tree = tree.root_at_midpoint()

    with open(args.output_pre +"_HOR_flank_dist_weighted.nwk", "w") as file:
        print(tree,file=file)

if __name__ == '__main__':
    main()
