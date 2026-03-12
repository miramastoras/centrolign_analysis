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

v3: Added --indel_window parameter. SNVs within --indel_window reference bases
of any indel in the pairwise cigar are excluded from Ti/Tv calculations.
Total nucleotide differences and nucleotide diversity remain unfiltered.
"""

import sys
import os
import re
import argparse
import numpy as np


def parse_cigar(cigar):
    parsed = []
    for m in re.finditer(r"(\d+)([HSNMIDX=])", cigar):
        parsed.append((m.group(2), int(m.group(1))))
    return parsed


def to_simple_cigar(cigar):
    simple = []
    for m in re.finditer(r"(\d+)([HSNMIDX=])", cigar):
        op = m.group(2)
        l = int(m.group(1))
        if op == "=" or op == "X":
            op = "M"
        if len(simple) != 0 and op == simple[-1][1]:
            simple[-1][0] += l
        else:
            simple.append([l, op])
    return "".join(str(l) + op for l, op in simple)


def count_aligned(cigar_ops):
    aligned = 0
    for op, l in cigar_ops:
        if op in ("M", "X", "="):
            aligned += l
    return aligned


def get_near_indel_flags(cigar_ops, window):
    """
    Walk the cigar and return a boolean numpy array of length equal to the
    number of X (mismatch) operations. True means that mismatch falls within
    `window` reference bases of an indel and should be excluded from Ti/Tv.

    Reference position advances for =, M, X, and D operations.
    Insertions (I) are anchored at the current reference position (between
    the preceding and following ref bases) and do not advance ref_pos.
    Deletions (D) are anchored at their start and end ref positions.
    """
    ref_pos = 0
    indel_anchors = []           # ref positions of indel events
    mismatch_ref_positions = []  # ref position of each X base

    for op, length in cigar_ops:
        if op in ("=", "M"):
            ref_pos += length
        elif op == "X":
            for i in range(length):
                mismatch_ref_positions.append(ref_pos + i)
            ref_pos += length
        elif op == "D":
            # Deletion in query: ref advances; anchor both edges
            indel_anchors.append(ref_pos)
            indel_anchors.append(ref_pos + length - 1)
            ref_pos += length
        elif op == "I":
            # Insertion in query: ref does not advance; anchor current ref_pos
            indel_anchors.append(ref_pos)
        # S, H, N: not modelled as indels here

    if not mismatch_ref_positions:
        return np.array([], dtype=bool)
    if not indel_anchors:
        return np.zeros(len(mismatch_ref_positions), dtype=bool)

    m_arr = np.array(mismatch_ref_positions, dtype=np.int64)
    i_arr = np.array(indel_anchors, dtype=np.int64)

    # For each mismatch, is there any indel anchor within `window` ref bases?
    near_indel = np.any(np.abs(m_arr[:, None] - i_arr[None, :]) <= window, axis=1)
    return near_indel


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Estimate nucleotide diversity from a SNP matrix and pairwise alignments, "
                    "optionally restricted to sample pairs within a distance threshold. "
                    "SNVs within --indel_window ref bp of any indel are excluded from Ti/Tv.")
    parser.add_argument("snp_mat", help="SNP matrix TSV file (from make_snp_mat)")
    parser.add_argument("pairwise_aln_prefix", help="Prefix path for pairwise cigar files, "
                        "e.g. /path/to/induced_pairwise_cigars/pairwise_cigar_")
    parser.add_argument("distance_csv", help="CSV file with columns sample1,sample2,distance")
    parser.add_argument("--threshold", type=float, default=0.2,
                        help="Only include pairs with distance < threshold (default: 0.2)")
    parser.add_argument("--indel_window", type=int, default=20,
                        help="Exclude SNVs within this many ref bp of any indel from Ti/Tv "
                             "(default: 20). Set to 0 to disable filtering.")
    args = parser.parse_args()

    prefix = os.path.abspath(args.pairwise_aln_prefix)
    cigar_dir = os.path.dirname(prefix)
    cigar_file_prefix = os.path.basename(prefix)

    # Load pairs from distance CSV: track all pairs seen, and those below threshold
    allowed_pairs = set()
    all_csv_pairs = set()
    with open(args.distance_csv) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split(',')
            s1, s2, dist = parts[0], parts[1], parts[2]
            all_csv_pairs.add(frozenset([s1, s2]))
            if float(dist) < args.threshold:
                allowed_pairs.add(frozenset([s1, s2]))

    print("Pairs within distance threshold {}: {}".format(args.threshold, len(allowed_pairs)),
          file=sys.stderr)

    # Read SNP matrix into memory: {sample: [allele_or_None, ...]}
    sample_alleles = {}
    with open(args.snp_mat) as fh:
        next(fh)  # discard header
        for line in fh:
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

    # Single pass over pairs: accumulate total_aligned, total_diffs (unfiltered),
    # and Ti/Tv counts (filtered by indel_window using cigar).
    #
    # total_diffs / total_aligned = nucleotide diversity (unfiltered, as before).
    # Ti/Tv uses only pairs where a cigar file was found; SNVs within indel_window
    # ref bp of any indel in that pair's cigar are excluded.
    total_aligned = 0
    total_diffs = 0
    pairs_found = 0
    pairs_missing = 0
    num_transitions = 0
    num_transversions = 0
    num_snvs_filtered = 0
    # Per-position flags: was there >=1 allowed pair with a (filtered) Ti/Tv here?
    ti_pos_flags = np.zeros(num_positions, dtype=bool)
    tv_pos_flags = np.zeros(num_positions, dtype=bool)

    for pair in allowed_pairs:
        s1, s2 = tuple(pair)
        both_in_matrix = s1 in sample_idx and s2 in sample_idx

        # --- Cigar lookup ---
        cigar_ops = None
        for f1, f2 in [(s1, s2), (s2, s1)]:
            fname = os.path.join(cigar_dir, "{}{}_{}.txt".format(cigar_file_prefix, f1, f2))
            if os.path.exists(fname):
                with open(fname) as fh:
                    cigar_ops = parse_cigar(fh.read().strip())
                total_aligned += count_aligned(cigar_ops)
                pairs_found += 1
                break
        if cigar_ops is None:
            pairs_missing += 1

        if not both_in_matrix:
            continue

        row_i = mat[sample_idx[s1]]
        row_j = mat[sample_idx[s2]]

        # All differing positions for this pair (unfiltered -> nucleotide diversity)
        valid_diff = (row_i != 0) & (row_j != 0) & (row_i != row_j)
        diff_positions = np.where(valid_diff)[0]
        total_diffs += len(diff_positions)

        if len(diff_positions) == 0:
            continue

        # --- Indel filtering for Ti/Tv ---
        if cigar_ops is not None and args.indel_window > 0:
            near_indel = get_near_indel_flags(cigar_ops, args.indel_window)
            if len(near_indel) != len(diff_positions):
                print(
                    "WARNING: cigar X count ({}) != SNP matrix differences ({}) "
                    "for pair {}, {}. Skipping indel filter for this pair.".format(
                        len(near_indel), len(diff_positions), s1, s2),
                    file=sys.stderr)
                keep_positions = diff_positions
            else:
                num_snvs_filtered += int(near_indel.sum())
                keep_positions = diff_positions[~near_indel]
        else:
            # No cigar found or window=0: include all diffs without indel filtering
            keep_positions = diff_positions

        if len(keep_positions) == 0:
            continue

        kept_a1 = row_i[keep_positions]
        kept_a2 = row_j[keep_positions]
        is_ti = TI_MATRIX[kept_a1, kept_a2]

        num_transitions += int(is_ti.sum())
        num_transversions += int((~is_ti).sum())
        ti_pos_flags[keep_positions[is_ti]] = True
        tv_pos_flags[keep_positions[~is_ti]] = True

    num_transition_events = int(ti_pos_flags.sum())
    num_transversion_events = int(tv_pos_flags.sum())

    print("Cigar files found: {}, missing: {}".format(pairs_found, pairs_missing),
          file=sys.stderr)
    print("SNVs filtered (within {} ref bp of indel): {}".format(
          args.indel_window, num_snvs_filtered), file=sys.stderr)
    print("total aligned pairs: {}".format(total_aligned))
    print("total nucleotide differences: {}".format(total_diffs))
    print("nucleotide diversity: {}".format(
          float(total_diffs) / total_aligned if total_aligned else 0))
    if num_transversions > 0:
        print("Ti/Tv ratio (differences): {}".format(
              float(num_transitions) / num_transversions))
    if num_transversion_events > 0:
        print("Ti/Tv ratio (events): {}".format(
              float(num_transition_events) / num_transversion_events))
