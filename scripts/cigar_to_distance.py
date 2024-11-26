#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 18:01:16 2023

Process pairwise alignments into a distance matrix and estimate a tree
using the neighbor joining algorithm

@author: Jordan
"""

import sys
import os
import re

def parse_cigar(cigar):
    parsed = []
    for m in re.finditer("(\d+)([HSMIDX=])", cigar):
        parsed.append((m.group(2), int(m.group(1))))
    return parsed

def cigar_to_dist(cigar, min_scale):

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

    if min_scale:
        return 1.0 - matches / min(ref_len, query_len)
    else:
        return 1.0 - (2.0 * matches) / (ref_len + query_len)

if __name__ == "__main__":

    if len(sys.argv) not in [2, 3]:
        print("usage: ./cigar_to_distance.py aln_dir [use_min_scale]", file = sys.stderr)
        exit(1)

    aln_dir = sys.argv[1]

    use_min_scale = False
    if len(sys.argv) == 3:
        use_min_scale = bool(int(sys.argv[2]))

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
        dist = cigar_to_dist(cigar, use_min_scale)
        with open(aln_dir+"pairwise_distance.csv", 'a') as f:
            print(samp1,samp2,dist,sep=",",file=f)
