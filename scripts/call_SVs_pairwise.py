#!/usr/bin/env python3

"""
Call SVs from centrolign pairwise cigar strings.
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
        prog='call_SVs_pairwise.py',
        description="""Call SVs from centrolign pairwise cigar strings. Returns bedPE file of SVs.""")

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

def cigar_to_sv_positions(cigar_ops, bedfile_prefix,len_diff_threshold=0.1, min_len=50,ref_name="ref",query_name="query"):
    """
    Iterate through CIGAR operations and return (ref, query) position pairs
    for insertions/deletions > 50 bp (or meeting adjacency rules).

    Adjacency Rules:
    - Include SV if I/D > 50 bp AND flanked by M on both sides.
    - Include SV if I/D adjacent to D/I AND (EITHER adjacent op is > 50bp) AND lengths differ by > len_diff_threshold (default 10%).
    """
    bedfile=bedfile_prefix+ref_name+"_"+query_name+".bed"
    #positions = []
    ref_pos = 0
    query_pos = 0

    def len_diff_above_threshold(a, b, threshold=0.1):
        """Return True if lengths differ by more than the given fractional threshold."""
        if a == 0 or b == 0:
            return False
        return abs(a - b) / max(a, b) > threshold

    for i, (op, length) in enumerate(cigar_ops):
        prev_op = cigar_ops[i - 1][0] if i > 0 else None
        prev_len = cigar_ops[i - 1][1] if i > 0 else None
        next_op = cigar_ops[i + 1][0] if i < len(cigar_ops) - 1 else None
        next_len = cigar_ops[i + 1][1] if i < len(cigar_ops) - 1 else None

        if op in "MX=":
            ref_pos += length
            query_pos += length

        elif op in "IHS":
            include = False
            # Case 1: flanked by matches on both sides and > 50 bp
            if length > min_len and prev_op == 'M' and next_op == 'M':
                include = True
            # Case 2: adjacent to D, EITHER is > min_len, and length diff > threshold (determines if true SV or just unaligned sequence)
            elif next_op == 'D' and (length > min_len or next_len > min_len) and \
                 len_diff_above_threshold(length, next_len, len_diff_threshold):
                include = True
            elif prev_op == 'D' and (length > min_len or prev_len > min_len) and \
                 len_diff_above_threshold(length, prev_len, len_diff_threshold):
                include = True

            if include:
                #positions.append((ref_name, ref_pos, ref_pos+1,  # end = ref_pos (exclusive)
                           #query_name, query_pos, query_pos + length, 'I'))  # end exclusive
                #print(ref_name, ref_pos, ref_pos + 1, query_name, query_pos, query_pos + length, "I",sep="\t")
                with open(bedfile,"a") as f:
                    print(ref_name, ref_pos, ref_pos + 1, query_name, query_pos, query_pos + length, "I", sep="\t",file=f)

            query_pos += length

        elif op == 'D':
            include = False
            # Case 1: flanked by matches on both sides and > 50 bp
            if length > min_len and prev_op == 'M' and next_op == 'M':
                include = True
            # Case 2: adjacent to I, EITHER > min_len, and large diff
            elif next_op == 'I' and (length > min_len or next_len > min_len) and \
                 len_diff_above_threshold(length, next_len, len_diff_threshold):
                include = True
            elif prev_op == 'I' and (length > min_len or prev_len > min_len) and \
                 len_diff_above_threshold(length, prev_len, len_diff_threshold):
                include = True

            if include:
                #positions.append((ref_name, ref_pos, ref_pos + length,  # end exclusive
                               #query_name, query_pos, query_pos+1, 'D'))  # end exclusive
                # print(ref_name, ref_pos, ref_pos + length, query_name, query_pos, query_pos + 1, "D", sep="\t")
                with open(bedfile, "a") as f:
                    print(ref_name, ref_pos, ref_pos + length, query_name, query_pos, query_pos + 1, "D", sep="\t",file=f)


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

    # for each cigar string, call SVs and write to bedfile
    for pair in sorted(filtered_pairs, key=lambda k: sorted(k)):
        # parse cigar string
        cigar = parse_cigar(open(filtered_pairs[pair]).read())

        print(f"\nProcessing: {pair[0]} vs {pair[1]} ({os.path.basename(filtered_pairs[pair])})",file=sys.stderr)
        ref=pair[0]
        query=pair[1]

        # call SVs
        cigar_to_sv_positions(cigar,args.bed_prefix,ref_name=ref,query_name=query)

    # skipping write to bed file for now

if __name__ == "__main__":
    main()
