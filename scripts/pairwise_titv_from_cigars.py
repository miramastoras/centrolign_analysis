#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute pairwise Ti/Tv ratio directly from pairwise CIGAR strings and per-sample
FASTA files, without requiring a SNP matrix.

For each pairwise_cigar_{s1}_{s2}.txt in the cigar directory:
  - Load s1 (reference) and s2 (query) FASTA sequences
  - Walk the alignment: at each X (mismatch) position, read the actual bases
  - Classify each mismatch as transition (A<->G, C<->T) or transversion
  - Optionally filter SNVs within --indel_window ref bp of any indel

FASTA files are expected at: {fasta_dir}/{sample}_{chr}_hor_array.fasta
Cigar files are expected as: {cigar_prefix}{s1}_{s2}.txt
  e.g. /path/to/induced_pairwise_cigars/pairwise_cigar_
  (s1 = reference, advances on =, X, D; s2 = query, advances on =, X, I)

Optionally, a distance CSV (sample1,sample2,distance) and --threshold can be
provided to restrict processing to pairs below the threshold.

Per-pair results are written as TSV to stdout.
Summary statistics are written to stderr.
"""

import sys
import os
import re
import argparse

TRANSITIONS = {frozenset(["A", "G"]), frozenset(["C", "T"])}


def parse_cigar(cigar):
    parsed = []
    for m in re.finditer(r"(\d+)([HSNMIDX=])", cigar):
        parsed.append((m.group(2), int(m.group(1))))
    return parsed


def load_fasta(path):
    """Load the first sequence from a FASTA file as an uppercase string."""
    seq = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if seq:
                    break
            else:
                seq.append(line.strip())
    return "".join(seq).upper()


def get_mismatch_bases(cigar_ops, seq1, seq2):
    """
    Walk the alignment and return a list of (base1, base2) for every X position.
    seq1 = reference (s1), seq2 = query (s2).
    """
    ref_pos = 0
    qry_pos = 0
    mismatches = []
    for op, length in cigar_ops:
        if op in ("=", "M"):
            ref_pos += length
            qry_pos += length
        elif op == "X":
            for i in range(length):
                b1 = seq1[ref_pos + i] if ref_pos + i < len(seq1) else "N"
                b2 = seq2[qry_pos + i] if qry_pos + i < len(seq2) else "N"
                mismatches.append((b1, b2))
            ref_pos += length
            qry_pos += length
        elif op == "D":
            ref_pos += length
        elif op == "I":
            qry_pos += length
    return mismatches


def get_near_indel_flags(cigar_ops, window):
    """
    Return a list of bools, one per X base, indicating whether that mismatch
    falls within `window` reference bases of any indel.

    Reference position advances for =, M, X, and D operations.
    Insertions (I) are anchored at the current reference position.
    Deletions (D) are anchored at their start and end ref positions.
    """
    ref_pos = 0
    indel_anchors = []
    mismatch_ref_positions = []

    for op, length in cigar_ops:
        if op in ("=", "M"):
            ref_pos += length
        elif op == "X":
            for i in range(length):
                mismatch_ref_positions.append(ref_pos + i)
            ref_pos += length
        elif op == "D":
            indel_anchors.append(ref_pos)
            indel_anchors.append(ref_pos + length - 1)
            ref_pos += length
        elif op == "I":
            indel_anchors.append(ref_pos)

    if not mismatch_ref_positions:
        return []
    if not indel_anchors:
        return [False] * len(mismatch_ref_positions)

    flags = []
    for mp in mismatch_ref_positions:
        near = any(abs(mp - anchor) <= window for anchor in indel_anchors)
        flags.append(near)
    return flags


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Compute pairwise Ti/Tv from CIGAR strings + FASTA sequences.")
    parser.add_argument("cigar_prefix",
                        help="Prefix path for cigar files, e.g. "
                             "/path/to/induced_pairwise_cigars/pairwise_cigar_")
    parser.add_argument("fasta_dir",
                        help="Directory containing {sample}_{chr}_hor_array.fasta files")
    parser.add_argument("chr",
                        help="Chromosome name used in FASTA filenames (e.g. chr14)")
    parser.add_argument("--indel_window", type=int, default=20,
                        help="Exclude SNVs within this many ref bp of any indel from Ti/Tv "
                             "(default: 20). Set to 0 to disable filtering.")
    parser.add_argument("--distance_csv",
                        help="Optional CSV (sample1,sample2,distance) to restrict pairs")
    parser.add_argument("--threshold", type=float, default=0.2,
                        help="Only process pairs with distance < threshold "
                             "(requires --distance_csv, default: 0.2)")
    args = parser.parse_args()

    # Build allowed-pair filter from distance CSV if provided
    allowed_pairs = None
    if args.distance_csv:
        allowed_pairs = set()
        with open(args.distance_csv) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split(",")
                if len(parts) < 3:
                    continue
                s1, s2, dist = parts[0], parts[1], parts[2]
                if float(dist) < args.threshold:
                    allowed_pairs.add(frozenset([s1, s2]))
        print("Pairs within distance threshold {}: {}".format(
              args.threshold, len(allowed_pairs)), file=sys.stderr)

    # Discover cigar files
    prefix = os.path.abspath(args.cigar_prefix)
    cigar_dir = os.path.dirname(prefix)
    cigar_file_prefix = os.path.basename(prefix)
    cigar_files = sorted(
        f for f in os.listdir(cigar_dir) if f.startswith(cigar_file_prefix)
    )

    fasta_cache = {}
    seen_pairs = set()  # frozensets of (s1, s2) already processed

    total_ti = 0
    total_tv = 0
    total_filtered = 0
    pairs_processed = 0
    pairs_skipped_threshold = 0
    pairs_missing_fasta = 0
    pairs_skipped_duplicate = 0

    # Per-pair TSV header
    print("sample1\tsample2\ttransitions\ttransversions\tti_tv_ratio\tsnvs_filtered")

    for fname in cigar_files:
        # Parse s1, s2: strip the file prefix and .txt suffix to get "{s1}_{s2}"
        name = fname[len(cigar_file_prefix):-len(".txt")]
        parts = name.split("_", 1)
        if len(parts) != 2:
            print("WARNING: could not parse sample names from {}".format(fname),
                  file=sys.stderr)
            continue
        s1, s2 = parts

        # Skip if the reverse-order file was already processed
        pair_key = frozenset([s1, s2])
        if pair_key in seen_pairs:
            pairs_skipped_duplicate += 1
            print("WARNING: duplicate cigar pair skipped (already seen): {}".format(fname),
                  file=sys.stderr)
            continue
        seen_pairs.add(pair_key)

        # Apply distance threshold filter
        if allowed_pairs is not None:
            if frozenset([s1, s2]) not in allowed_pairs:
                pairs_skipped_threshold += 1
                continue

        # Load FASTAs (cached)
        fasta_path1 = os.path.join(args.fasta_dir, "{}_{}_hor_array.fasta".format(s1, args.chr))
        fasta_path2 = os.path.join(args.fasta_dir, "{}_{}_hor_array.fasta".format(s2, args.chr))

        for path in (fasta_path1, fasta_path2):
            if path not in fasta_cache:
                if os.path.exists(path):
                    fasta_cache[path] = load_fasta(path)
                else:
                    fasta_cache[path] = None

        seq1 = fasta_cache[fasta_path1]
        seq2 = fasta_cache[fasta_path2]

        if seq1 is None or seq2 is None:
            missing = fasta_path1 if seq1 is None else fasta_path2
            print("WARNING: FASTA not found: {}".format(missing), file=sys.stderr)
            pairs_missing_fasta += 1
            continue

        # Parse cigar
        cigar_path = os.path.join(cigar_dir, fname)
        with open(cigar_path) as fh:
            cigar_str = fh.read().strip()
        cigar_ops = parse_cigar(cigar_str)

        # Get mismatch bases and indel-proximity flags
        mismatches = get_mismatch_bases(cigar_ops, seq1, seq2)
        if not mismatches:
            print("{}\t{}\t0\t0\tNA\t0".format(s1, s2))
            pairs_processed += 1
            continue

        if args.indel_window > 0:
            near_indel = get_near_indel_flags(cigar_ops, args.indel_window)
        else:
            near_indel = [False] * len(mismatches)

        ti = 0
        tv = 0
        filtered = 0
        for i, (b1, b2) in enumerate(mismatches):
            if near_indel[i]:
                filtered += 1
                continue
            if b1 == "N" or b2 == "N":
                continue
            if frozenset([b1, b2]) in TRANSITIONS:
                ti += 1
            else:
                tv += 1

        ratio = "{:.6f}".format(ti / tv) if tv > 0 else "NA"
        print("{}\t{}\t{}\t{}\t{}\t{}".format(s1, s2, ti, tv, ratio, filtered))

        total_ti += ti
        total_tv += tv
        total_filtered += filtered
        pairs_processed += 1

    # Summary
    overall_ratio = "{:.6f}".format(total_ti / total_tv) if total_tv > 0 else "NA"
    print("Pairs processed: {}".format(pairs_processed), file=sys.stderr)
    if pairs_skipped_duplicate > 0:
        print("Pairs skipped (duplicate forward/reverse): {}".format(pairs_skipped_duplicate),
              file=sys.stderr)
    if allowed_pairs is not None:
        print("Pairs skipped (above threshold): {}".format(pairs_skipped_threshold),
              file=sys.stderr)
    print("Pairs missing FASTA: {}".format(pairs_missing_fasta), file=sys.stderr)
    print("Total transitions: {}".format(total_ti), file=sys.stderr)
    print("Total transversions: {}".format(total_tv), file=sys.stderr)
    print("Overall Ti/Tv ratio: {}".format(overall_ratio), file=sys.stderr)
    print("SNVs filtered (within {} ref bp of indel): {}".format(
          args.indel_window, total_filtered), file=sys.stderr)
