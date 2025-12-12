#!/usr/bin/env python3

"""
Call short indels from centrolign pairwise cigar strings.
"""

import sys
import os
import re
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from itertools import combinations


def arg_parser():
    '''
    Parses command line arguments with argparse
    '''
    parser = argparse.ArgumentParser(
        prog='call_short_indels_pairwise.py',
        description="""Call short indels from centrolign pairwise cigar strings. Returns bedPE file of short indels.""")

    parser.add_argument("-c", "--cigars",
                        required=True,
                        help="Folder containing cigar strings. Must be named as {args.cigars}{smp1}_{smp2}.txt")
    parser.add_argument("-s", "--samples",
                        required=True,
                        help="Txt file containing list of samples to run script on.")
    parser.add_argument("-o", "--bed_prefix",
                        required=True,
                        help="Prefix / filepath to write output bed files to")

    return parser.parse_args()

def map_to_samples(fps):
    suff_regex = "([a-zA-Z0-9]+.[0-9])_([a-zA-Z0-9]+.[0-9]).txt$"

    sample_map = {}
    for fp in fps:
        m = re.search(suff_regex, fp)
        assert (m is not None)
        key = tuple([m.group(1), m.group(2)])
        sample_map[key] = fp

    return sample_map

def cigar_to_indel_positions(cigar_ops, bedfile_prefix, max_len=49,ref_name="ref",query_name="query"):
    """
    Iterate through CIGAR operations and return (ref, query) position pairs
    for insertions/deletions < 50 bp
    """
    bedfile=bedfile_prefix+ref_name+"_"+query_name+".bed"
    #positions = []
    ref_pos = 0
    query_pos = 0

    def percent_diff(a, b):
        """Calculate percent difference between two values"""
        if a == 0 or b == 0:
            return False
        return abs(a - b) / max(a, b)

    for i, (op, length) in enumerate(cigar_ops):
        # if prev_op is the beginning or next_op is the end, set it to "M"
        prev_op = cigar_ops[i - 1][0] if i > 0 else "M"
        prev_len = cigar_ops[i - 1][1] if i > 0 else None
        next_op = cigar_ops[i + 1][0] if i < len(cigar_ops) - 1 else "M"
        next_len = cigar_ops[i + 1][1] if i < len(cigar_ops) - 1 else None

        # Don't write SNPs but iterate the position counters
        if op in "MX=":
            ref_pos += length
            query_pos += length

        elif op in "IHS":

            if length <= max_len: # only write short indels < 50 bp
                adj_id_diff = ""

                # Case 1: flanked by matches on both sides
                if prev_op in "MX=" and next_op in "MX=":
                    adj_id_diff = -1

                # Case 2: I adjacent to a D
                elif next_op == 'D' :
                    adj_id_diff = percent_diff(length,next_len)

                elif prev_op == 'D' :
                    adj_id_diff = percent_diff(length,prev_len)

                with open(bedfile,"a") as f:
                    print(ref_name, ref_pos, ref_pos + 1, query_name, query_pos, query_pos + length, "I", adj_id_diff, sep="\t",file=f)

            query_pos += length

        elif op == 'D':

            if length <= max_len:

                adj_id_diff=""
                # Case 1: flanked by matches on both sides and >= 50 bp
                if prev_op in "MX=" and next_op in "MX=":
                    adj_id_diff = -1

                # Case 2: adjacent to I, EITHER is <= max_len, and large diff
                elif next_op in 'IHS' :
                    adj_id_diff = percent_diff(length,next_len)

                elif prev_op in 'IHS' :
                    adj_id_diff = percent_diff(length,prev_len)

                with open(bedfile, "a") as f:
                    print(ref_name, ref_pos, ref_pos + length, query_name, query_pos, query_pos + 1, "D", adj_id_diff, sep="\t",file=f)

            ref_pos += length

        else:
            assert (False)

    return

def parse_cigar(cigar):
    '''
    :param cigar: string
    :return: ordered list of cigar operations: [(op,length)]
    '''
    parsed = []
    for m in re.finditer("(\d+)([HSMIDX=])", cigar):
        parsed.append((m.group(2), int(m.group(1))))
    return parsed

# on real data double check there aren't any DID or IDI

def main():
    # parse command line arguments
    args = arg_parser()

    # read in list of samples
    with open(args.samples) as f:
        samples = [line.strip() for line in f if line.strip()]
    print(f"Loaded {len(samples)} samples from {args.samples}",file=sys.stderr)

    # generate all pairwise combinations (unordered)
    valid_pairs = set(tuple(sorted(pair)) for pair in combinations(samples, 2))

    # read in all cigar strings from input directory
    cigar_prefix = os.path.abspath(args.cigars)
    cigar_dir = os.path.dirname(cigar_prefix)
    cigar_file_prefix = os.path.basename(cigar_prefix)

    cigar_files = [os.path.join(cigar_dir, f) for f in os.listdir(cigar_dir) if f.startswith(cigar_file_prefix)]

    # map files to sample pairs
    samples_to_induced = map_to_samples(cigar_files)

    # filter only those pairs present in the provided sample list
    filtered_pairs = {pair: fp for pair, fp in samples_to_induced.items() if tuple(sorted(pair)) in valid_pairs}

    print(f"Running on {len(filtered_pairs)} pairwise combinations:",file=sys.stderr)

    # for each cigar string, call short indels and write to bedfile
    for pair in sorted(filtered_pairs, key=lambda k: sorted(k)):
        # parse cigar string
        cigar = parse_cigar(open(filtered_pairs[pair]).read())

        print(f"\nProcessing: {pair[0]} vs {pair[1]} ({os.path.basename(filtered_pairs[pair])})",file=sys.stderr)
        ref=pair[0]
        query=pair[1]

        # call short indels
        cigar_to_indel_positions(cigar,args.bed_prefix,ref_name=ref,query_name=query)

    # skipping write to bed file for now

if __name__ == "__main__":
    main()
