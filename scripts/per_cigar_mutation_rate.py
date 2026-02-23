#!/usr/bin/env python3

"""
Calculate mutation rates per CIGAR string from centrolign pairwise alignments.
Combines SNV, short indel, and SV counting in a single pass through each CIGAR.

Output columns per pair:
  - n_svs, n_svs_per_avg_len
  - n_snvs, n_snvs_per_aligned_base, n_snvs_per_aligned_base_per_avg_len
  - n_short_indels, n_short_indels_per_aligned_base, n_short_indels_per_aligned_base_per_avg_len

For short indels (I/D <= 49bp) and SVs (I/D > 49bp), adjacent I/D ops are
only counted if adj_id_diff > 0.1. Ops flanked by matches on both sides are
always counted.
"""

import sys
import os
import re
import argparse
import pandas as pd
from pathlib import Path



def arg_parser():
    parser = argparse.ArgumentParser(
        prog='per_cigar_mutation_rate.py',
        description="Calculate mutation rates per CIGAR string from centrolign pairwise alignments.")

    parser.add_argument("-c", "--contig_map_csv",
                        required=True,
                        help="CSV file with columns: sample1, sample1_contig, sample2, sample2_contig, pairwise_cigar_path")
    parser.add_argument("-o", "--output",
                        required=True,
                        help="Output TSV file path")
    parser.add_argument("-a", "--asat_bed_dir",
                        required=True,
                        help="Directory with per-sample ASAT bed files (*_asat_arrays.bed)")
    parser.add_argument("-chr", "--chrom",
                        required=True,
                        help="Chromosome name (e.g. chr1)")

    return parser.parse_args()


def parse_cigar(cigar):
    '''
    :param cigar: string
    :return: ordered list of cigar operations: [(op, length)]
    '''
    parsed = []
    for m in re.finditer(r"(\d+)([HSMIDX=])", cigar):
        parsed.append((m.group(2), int(m.group(1))))
    return parsed



def percent_diff(a, b):
    """Calculate percent difference between two values."""
    if a == 0 or b == 0:
        return False
    return abs(a - b) / max(a, b)


def passes_adj_filter(op, length, prev_op, prev_len, next_op, next_len):
    """
    Return True if this I/D op should be counted based on adjacency rules.

    Case 1 (flanked by matches on both sides): always count.
    Case 2 (adjacent to opposite indel type): only count if adj_id_diff > 0.1.
    Otherwise: do not count.
    """
    # Case 1: flanked by matches on both sides
    if prev_op in "MX=" and next_op in "MX=":
        return True

    # Case 2: adjacent to opposite indel type
    if op in "IHS":
        if next_op == "D":
            pd_val = percent_diff(length, next_len)
            return pd_val is not False and pd_val > 0.1
        elif prev_op == "D":
            pd_val = percent_diff(length, prev_len)
            return pd_val is not False and pd_val > 0.1

    elif op == "D":
        if next_op in "IHS":
            pd_val = percent_diff(length, next_len)
            return pd_val is not False and pd_val > 0.1
        elif prev_op in "IHS":
            pd_val = percent_diff(length, prev_len)
            return pd_val is not False and pd_val > 0.1

    return False


def compute_mutation_counts(cigar_ops):
    """
    Single pass through CIGAR ops to count SNVs, short indels, SVs, and aligned bases.

    Returns: (n_snvs, n_short_indels, n_svs, aligned_bases)
    """
    n_snvs = 0
    aligned_bases = 0
    n_short_indels = 0
    n_svs = 0

    ref_pos = 0
    query_pos = 0

    for i, (op, length) in enumerate(cigar_ops):
        prev_op = cigar_ops[i - 1][0] if i > 0 else "M"
        prev_len = cigar_ops[i - 1][1] if i > 0 else None
        next_op = cigar_ops[i + 1][0] if i < len(cigar_ops) - 1 else "M"
        next_len = cigar_ops[i + 1][1] if i < len(cigar_ops) - 1 else None

        if op in "MX=":
            if op == "X":
                n_snvs += length
            aligned_bases += length
            ref_pos += length
            query_pos += length

        elif op in "IHS":
            if op == "I":
                if length <= 49:
                    if passes_adj_filter(op, length, prev_op, prev_len, next_op, next_len):
                        n_short_indels += length
                else:  # length > 49
                    if passes_adj_filter(op, length, prev_op, prev_len, next_op, next_len):
                        n_svs += length
            query_pos += length

        elif op == "D":
            if length <= 49:
                if passes_adj_filter(op, length, prev_op, prev_len, next_op, next_len):
                    n_short_indels += length
            else:  # length > 49
                if passes_adj_filter(op, length, prev_op, prev_len, next_op, next_len):
                    n_svs += length
            ref_pos += length

        else:
            assert False, f"Unknown CIGAR op: {op}"

    return n_snvs, n_short_indels, n_svs, aligned_bases


def load_asat_bounds(asat_bed_dir):
    """Load ASAT bed files and return dict of (sample, chr) -> (start, end)."""
    asat_dir = Path(asat_bed_dir)
    dfs = []
    for bed_file in asat_dir.glob("*_asat_arrays.bed"):
        sample = bed_file.name.replace("_asat_arrays.bed", "")
        df = pd.read_csv(bed_file, sep="\t", header=None,
                         names=["contig", "start", "end", "chr"])
        df["sample"] = sample
        dfs.append(df)
    asat_df = pd.concat(dfs, ignore_index=True)

    asat_bounds = {}
    for row in asat_df.itertuples(index=False):
        contig, start, end, chr_, sample = row
        asat_bounds[(sample, chr_)] = (start, end)

    return asat_bounds


def main():
    args = arg_parser()

    # Load ASAT bounds for array lengths
    print("Loading ASAT beds...", file=sys.stderr)
    asat_bounds = load_asat_bounds(args.asat_bed_dir)
    print(f"Loaded ASAT bounds for {len(asat_bounds)} (sample, chr) pairs", file=sys.stderr)

    # Read contig map CSV
    contig_map = pd.read_csv(args.contig_map_csv)
    print(f"Loaded {len(contig_map)} pairs from {args.contig_map_csv}", file=sys.stderr)

    # Process each pair
    results = []
    for _, row in contig_map.iterrows():
        sample1 = row["sample1"]
        sample2 = row["sample2"]
        cigar_path = row["pairwise_cigar_path"]

        if not os.path.exists(cigar_path):
            print(f"  WARNING: CIGAR file not found: {cigar_path}, skipping", file=sys.stderr)
            continue

        print(f"  Processing: {sample1} vs {sample2}", file=sys.stderr)

        cigar_ops = parse_cigar(open(cigar_path).read())
        n_snvs, n_short_indels, n_svs, aligned_bases = compute_mutation_counts(cigar_ops)

        # Get array lengths for the two samples
        s1_bounds = asat_bounds.get((sample1, args.chrom))
        s2_bounds = asat_bounds.get((sample2, args.chrom))

        if s1_bounds is None or s2_bounds is None:
            print(f"    WARNING: missing ASAT bounds for {sample1} or {sample2} "
                  f"on {args.chrom}, skipping", file=sys.stderr)
            continue

        s1_array_len = s1_bounds[1] - s1_bounds[0]
        s2_array_len = s2_bounds[1] - s2_bounds[0]
        avg_array_len = (s1_array_len + s2_array_len) / 2

        out_row = {
            "sample1": sample1,
            "sample2": sample2,
            "chr": args.chrom,
            "avg_array_len": avg_array_len,
            "aligned_bases": aligned_bases,
            "n_snvs": n_snvs,
            "n_snvs_per_aligned_base": (n_snvs / aligned_bases
                                        if aligned_bases > 0 else 0),
            "n_snvs_per_aligned_base_per_avg_len": (n_snvs / aligned_bases / avg_array_len
                                                    if aligned_bases > 0 and avg_array_len > 0
                                                    else 0),
            "n_short_indels_bases": n_short_indels,
            "n_short_indels_per_aligned_base": (n_short_indels / aligned_bases
                                                if aligned_bases > 0 else 0),
            "n_short_indels_per_aligned_base_per_avg_len": (n_short_indels / aligned_bases / avg_array_len
                                                            if aligned_bases > 0 and avg_array_len > 0
                                                            else 0),
            "n_sv_bases": n_svs,
            "n_svs_per_avg_len": (n_svs / avg_array_len
                                  if avg_array_len > 0 else 0),
        }
        results.append(out_row)

    # Write output
    df = pd.DataFrame(results)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"\nWrote {len(df)} rows to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
