#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 18:01:16 2023

Process pairwise alignments into a distance matrix

@author: Jordan
"""

import sys
import os
import re
import argparse

def parse_cigar(cigar):
    parsed = []
    for m in re.finditer("(\d+)([HSMIDX=])", cigar):
        parsed.append((m.group(2), int(m.group(1))))
    return parsed


def cigar_to_dist_method3(cigar):
    '''
    dist = mismatches / (matches+mismatches)
    # this formula doesn't work for arrays with a low # of aligned bases ie 210=2359579D3069392I168=
    '''
    matches=0
    mismatches=0

    for op, op_len in cigar:
        print(op,op_len)
        if op == "X":
            mismatches += op_len
        elif op in "M=" :
            matches+= op_len
        else:
            mismatches += 1
    print(mismatches, matches)
    return mismatches / (matches + mismatches)

def cigar_to_dist_method2(cigar):
    '''
    Does not accept cigar strings that use M instead of X or =
    Calculates distance as follows:
    Distance = 1 - (proportion aligned) * (matches / (matches + mismatches))
    Proportion aligned = (ref aligned + query aligned / ref total + query total)
    Proportion aligned = ((matches+mismatches)*2)  / (((matches+mismatches)*2) + insertions + deletions)
    '''
    print("calculating distance 2")
    matches = 0
    mismatches=0
    deletions=0
    insertions=0

    for op, op_len in cigar:
        if op == "X":
            mismatches += op_len
        elif op == "=":
            matches += op_len
        elif op == "D":
            deletions+= op_len
        elif op in "I":
            insertions += op_len
        else:
            assert (False)

    # if there are no matches, distance is 1 (avoid div by zero error)
    if matches == 0:
        return 1.0

    # Proportion aligned = (ref aligned + query aligned / ref total + query total)
    prop_aligned = ((matches+mismatches)*2) / (((matches+mismatches)*2) + insertions + deletions)
    print(matches,mismatches,insertions,deletions)
    return 1 - (prop_aligned * (matches / (matches + mismatches)))

def cigar_to_dist_method1(cigar, min_scale):
    '''
    1.0 - (2.0 * matches)
       ---------------
     (ref_len + query_len)
    '''
    print("calculating distance 1 - Jordan's original formula")
    query_len = 0
    ref_len = 0
    matches = 0
    for op, op_len in cigar:
        if op in "MX=":
            query_len += op_len
            ref_len += op_len
            if op != "X":
                matches += op_len
        elif op == "D":
            ref_len += op_len
        elif op in "IHS":
            query_len += op_len
        else:
            assert(False)
    #print(matches, ref_len,query_len)
    if min_scale:
        return 1.0 - matches / min(ref_len, query_len)
    else:
        return 1.0 - (2.0 * matches) / (ref_len + query_len)

def arg_parser():
    '''
    Parses command line arguments with argparse
    '''
    parser = argparse.ArgumentParser(
        prog='cigar_to_distance.py',
        description="""converts cigars to distances""")

    parser.add_argument("-a", "--aln_dir",
                        required=True,
                        help="Directory containing pairwise cigar strings. They must be named as pairwise_cigar_smp1_smp2.txt")
    parser.add_argument("-d", "--dist_metric",
                        required=True,
                        help="Calculate distance using formula 1, 2 or 3.")
    parser.add_argument("-o", "--output_prefix",
                        required=True,
                        help="Prefix to write csv output file")
    parser.add_argument("--use_min_scale",
                        action='store_true',
                        help="Use min scaling with formula 1")

    return parser.parse_args()

if __name__ == "__main__":

    # parse command line arguments
    args = arg_parser()

    aln_dir = args.aln_dir

    use_min_scale = False

    if args.use_min_scale:
        use_min_scale=True

    mat = {}

    fps = os.listdir(aln_dir)

    for i in range(len(fps)):
        if (i + 1) % 100 == 0:
            print("processed {} of {} files".format(i + 1, len(fps)), file = sys.stderr)
        fp = fps[i]
        with open(os.path.join(aln_dir, fp)) as f:
            cigar = parse_cigar(f.read().strip())

        m = re.search("([0-9a-zA-Z.]+)_([0-9a-zA-Z.]+)(_cigar)?.txt", fp)
        assert(m is not None)
        samp1 = m.group(1)
        samp2 = m.group(2)

        print(int(args.dist_metric))
        if int(args.dist_metric) == 1:
            dist = cigar_to_dist_method1(cigar, use_min_scale)
        elif int(args.dist_metric) == 2:
            dist = cigar_to_dist_method2(cigar)
        elif int(args.dist_metric) == 3:
            dist = cigar_to_dist_method3(cigar)
        elif int(args.dist_metric) == 4:
            dist = cigar_to_dist_method4(cigar)
        else:
            print("Invalid entry for --dist_method. Options are 1, 2, 3, 4.", file=sys.stderr)
            exit(1)
        dist = "{:.20f}".format(dist)
        with open(args.output_prefix+"pairwise_distance.csv", 'a') as f:
            print(samp1,samp2,dist,sep=",",file=f)

#(( 2 x matches) / ( (2 x (matches + mismatches)) + insertions + deletions))

#2: 1 - ((( 2 x (matches+ mismatches)) / ( (2 x (matches + mismatches)) + insertions + deletions)) x (matches/(mismatches x matches)))
