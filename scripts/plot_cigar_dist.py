#!/usr/bin/env python3
"""
This script takes all cigar strings in an input directory, and plots the # mismatches vs proportion aligned
"""

import sys
import os
import re
import argparse
import matplotlib.pyplot as plt

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

    parser.add_argument("-l", "--label",
                        required=True,
                        help="Label for plot title")

    parser.add_argument("-o", "--output_prefix",
                        required=True,
                        help="Prefix to write csv output file")

    return parser.parse_args()

def calc_mismatch_prop_aligned(cigar):
    matches = 0
    mismatches = 0
    deletions = 0
    insertions = 0

    for op, op_len in cigar:
        if op == "X":
            mismatches += op_len
        elif op == "=":
            matches += op_len
        elif op == "D":
            deletions += op_len
        elif op in "I":
            insertions += op_len
        else:
            assert (False)

    # Proportion aligned = (ref aligned + query aligned / ref total + query total)
    prop_aligned = ((matches + mismatches) * 2) / (((matches + mismatches) * 2) + insertions + deletions)

    return mismatches , prop_aligned

def parse_cigar(cigar):
    parsed = []
    for m in re.finditer("(\d+)([HSMIDX=])", cigar):
        parsed.append((m.group(2), int(m.group(1))))
    return parsed


if __name__ == "__main__":

    # parse command line arguments
    args = arg_parser()

    aln_dir = args.aln_dir

    mat = {}

    fps = os.listdir(aln_dir)

    mismatches=[]
    prop_aln=[]

    for i in range(len(fps)):
        if (i + 1) % 100 == 0:
            print("processed {} of {} files".format(i + 1, len(fps)), file=sys.stderr)
        fp = fps[i]
        with open(os.path.join(aln_dir, fp)) as f:
            cigar = parse_cigar(f.read().strip())

        m , p = calc_mismatch_prop_aligned(cigar)

        mismatches.append(m)
        prop_aln.append(p)

    plt.scatter(prop_aln, mismatches, s=10,c='blue', alpha=0.5,edgecolors='none')
    plt.title(args.label)
    plt.xlabel('Proportion aligned')
    plt.ylabel('# of mismatches')

    plt.savefig(args.output_prefix + '.png', dpi=300, bbox_inches='tight')
