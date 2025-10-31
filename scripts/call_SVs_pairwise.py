#!/usr/bin/env python3

"""
Call SVs from centrolign pairwise cigar strings.
"""

import sys
import os
import re

def cigar_to_positions(cigar_ops, len_diff_threshold=0.1, min_len=50,ref_name="ref",query_name="query"):
    """
    Iterate through CIGAR operations and return (ref, query) position pairs
    for insertions/deletions > 50 bp (or meeting adjacency rules).

    Adjacency Rules:
    - Include SV if I/D > 50 bp AND flanked by M on both sides.
    - Include SV if I/D adjacent to D/I AND (EITHER adjacent op is > 50bp) AND lengths differ by > len_diff_threshold (default 10%).
    """

    positions = []
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

        if op == 'M':
            ref_pos += length
            query_pos += length

        elif op == 'I':
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
                positions.append((ref_name, ref_pos, ref_pos+1,  # end = ref_pos (exclusive)
                           query_name, query_pos, query_pos + length, 'I'))  # end exclusive

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
                positions.append((ref_name, ref_pos, ref_pos + length,  # end exclusive
                               query_name, query_pos, query_pos+1, 'D'))  # end exclusive
            ref_pos += length

        else:
            assert (False)

    return positions

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

    #cigar="15M400I7M150I156D8M"
    cigar="15M600D300I2M"
    cigar="150M60I70D100M"
    cigar="1500D2000I2M300D2M"
    # read in cigar from folder
    # write output bedfiles
    
    # plot length distribution
    #
    parsed=parse_cigar(cigar)
    #print(parsed)
    sv_positions = cigar_to_positions(parsed)
    print(sv_positions)

if __name__ == "__main__":
    main()
