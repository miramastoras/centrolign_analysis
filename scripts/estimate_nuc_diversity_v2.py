#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 13:53:20 2023

Use a SNP matrix from the make_snp_mat script and a set of pairwise
alignments to estimate quantities analogous to nucleotide diversity and
transition-transversion ratio

@author: Jordan

Updated on Wed Mar 11 2026

@author: Mira
"""

import sys
import os
import re
import argparse

def parse_cigar(cigar):
    parsed = []
    for m in re.finditer("(\d+)([HSNMIDX=])", cigar):
        parsed.append((m.group(2), int(m.group(1))))
    return parsed

def to_simple_cigar(cigar):
    simple = []
    for m in re.finditer("(\d+)([HSNMIDX=])", cigar):
        op = m.group(2)
        l = int(m.group(1))
        if op == "=" or op == "X":
            op = "M"
        if len(simple) != 0 and op == simple[-1][1]:
            simple[-1][0] += l
        else:
            simple.append([l, op])
    return "".join(str(l) + op for l, op in simple)


def count_aligned(cigar):
    aligned = 0
    for op, l in cigar:
        if op == "M" or op == "X" or op == "=":
            aligned += l

    return aligned

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Estimate nucleotide diversity from a SNP matrix and pairwise alignments, "
                    "optionally restricted to sample pairs within a distance threshold.")
    parser.add_argument("snp_mat", help="SNP matrix TSV file (from make_snp_mat)")
    parser.add_argument("pairwise_aln_prefix", help="Prefix path for pairwise cigar files, "
                        "e.g. /path/to/induced_pairwise_cigars/pairwise_cigar_")
    parser.add_argument("distance_csv", help="CSV file with columns sample1,sample2,distance")
    parser.add_argument("--threshold", type=float, default=0.2,
                        help="Only include pairs with distance < threshold (default: 0.2)")
    args = parser.parse_args()

    prefix = os.path.abspath(args.pairwise_aln_prefix)
    cigar_dir = os.path.dirname(prefix)
    cigar_file_prefix = os.path.basename(prefix)

    # Load pairs from distance CSV: track all pairs seen, and those below threshold
    allowed_pairs = set()
    all_csv_pairs = set()
    with open(args.distance_csv) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(',')
            s1, s2, dist = parts[0], parts[1], parts[2]
            all_csv_pairs.add(frozenset([s1, s2]))
            if float(dist) < args.threshold:
                allowed_pairs.add(frozenset([s1, s2]))

    print("Pairs within distance threshold {}: {}".format(args.threshold, len(allowed_pairs)), file=sys.stderr)

    # Sum aligned bases from cigar files for allowed pairs only
    total_aligned = 0
    pairs_found = 0
    pairs_missing = 0
    for pair in allowed_pairs:
        s1, s2 = tuple(pair)
        found = False
        for f1, f2 in [(s1, s2), (s2, s1)]:
            fname = os.path.join(cigar_dir, "{}{}_{}.txt".format(cigar_file_prefix, f1, f2))
            if os.path.exists(fname):
                with open(fname) as f:
                    cigar = parse_cigar(f.read().strip())
                    total_aligned += count_aligned(cigar)
                pairs_found += 1
                found = True
                break
        if not found:
            pairs_missing += 1

    print("Cigar files found: {}, missing: {}".format(pairs_found, pairs_missing), file=sys.stderr)

    # Read SNP matrix into memory: {sample: [allele_or_None, ...]}
    sample_alleles = {}
    with open(args.snp_mat) as f:
        # discard the header
        next(f)
        for line in f:
            tokens = line.strip().split()
            sample = tokens[0]
            sample_alleles[sample] = [None if t == "?" else t for t in tokens[1:]]

    num_positions = len(next(iter(sample_alleles.values()))) if sample_alleles else 0

    # Warn about matrix pairs absent from the CSV
    matrix_samples = sorted(sample_alleles.keys())
    missing_from_csv = []
    for i in range(len(matrix_samples)):
        for j in range(i + 1, len(matrix_samples)):
            pair = frozenset([matrix_samples[i], matrix_samples[j]])
            if pair not in all_csv_pairs:
                missing_from_csv.append((matrix_samples[i], matrix_samples[j]))
    if missing_from_csv:
        print("WARNING: {} pair(s) found in SNP matrix but absent from distance CSV:".format(
            len(missing_from_csv)), file=sys.stderr)
        for s1, s2 in missing_from_csv:
            print("  {} , {}".format(s1, s2), file=sys.stderr)

    import numpy as np

    # Encode alleles as integers; 0 = missing (?)
    ENC = {"A": 1, "C": 2, "G": 3, "T": 4}

    # Build encoded matrix for samples that appear in at least one allowed pair
    samples_in_pairs = sorted(
        set(s for pair in allowed_pairs for s in pair) & set(sample_alleles.keys())
    )
    sample_idx = {s: i for i, s in enumerate(samples_in_pairs)}

    mat = np.zeros((len(samples_in_pairs), num_positions), dtype=np.int8)
    for sample, idx in sample_idx.items():
        for pos, allele in enumerate(sample_alleles[sample]):
            if allele is not None:
                mat[idx, pos] = ENC.get(allele, 0)

    # Transition lookup: TI_MATRIX[enc_a1, enc_a2] == True iff (a1,a2) is a transition
    # Encoding: A=1 C=2 G=3 T=4
    TI_MATRIX = np.zeros((5, 5), dtype=bool)
    for a1, a2 in [(1, 3), (3, 1), (2, 4), (4, 2)]:  # A-G, G-A, C-T, T-C
        TI_MATRIX[a1, a2] = True

    num_transitions = 0
    num_transversions = 0
    # Per-position flags: was there ≥1 allowed pair with a Ti/Tv at this position?
    ti_pos_flags = np.zeros(num_positions, dtype=bool)
    tv_pos_flags = np.zeros(num_positions, dtype=bool)

    for pair in allowed_pairs:
        s1, s2 = tuple(pair)
        if s1 not in sample_idx or s2 not in sample_idx:
            continue
        row_i = mat[sample_idx[s1]]
        row_j = mat[sample_idx[s2]]

        # Positions where both samples have a non-missing allele and they differ
        valid_diff = (row_i != 0) & (row_j != 0) & (row_i != row_j)

        # Classify each differing position as Ti or Tv using fancy indexing
        is_ti = TI_MATRIX[row_i, row_j]  # shape (num_positions,)
        ti_mask = valid_diff & is_ti
        tv_mask = valid_diff & ~is_ti

        num_transitions += int(ti_mask.sum())
        num_transversions += int(tv_mask.sum())
        ti_pos_flags |= ti_mask
        tv_pos_flags |= tv_mask

    num_transition_events = int(ti_pos_flags.sum())
    num_transversion_events = int(tv_pos_flags.sum())

    total_diffs = num_transitions + num_transversions
    print("total aligned pairs: {}".format(total_aligned))
    print("total nucleotide differences: {}".format(total_diffs))
    print("nucleotide diversity: {}".format(float(total_diffs) / total_aligned if total_aligned else 0))
    if num_transversions > 0:
        print("Ti/Tv ratio (differences): {}".format(float(num_transitions) / num_transversions))
    if num_transversion_events > 0:
        print("Ti/Tv ratio (events): {}".format(float(num_transition_events) / num_transversion_events))
