#!/usr/bin/env python3
"""
Per-chromosome permutation test: does the top 10% local identity exclusive
region (top 10% minus CDR) have an elevated short indel rate compared to the
rest of the array (outside both CDR and top 10%)?

test_region      = top_10%_local_identity - CDR
comparison_region = complement(CDR ∪ top_10%_local_identity)

Parallelizes across chromosomes using multiprocessing.

Usage:
    python localID_minus_CDR_permutation_per_chrom_short_indels.py --n-reps 10000 --n-workers 8 --output-dir /path/to/output
"""

import argparse
import os
import re
import glob
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool

import numpy as np
import pandas as pd


# ── CIGAR helpers ────────────────────────────────────────────────────────────

def parse_cigar(cigar):
    return [
        (int(length), op)
        for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    ]

def reverse_cigar(parsed):
    converted = []
    for length, op in parsed:
        if op == "I":
            converted.append((length, "D"))
        elif op == "D":
            converted.append((length, "I"))
        else:
            converted.append((length, op))
    return converted

def build_aligned_intervals_from_ops(cigar_ops, start_pos):
    intervals = []
    pos = start_pos
    for length, op in cigar_ops:
        if op in ("M", "=", "X"):
            intervals.append((pos, pos + length))
            pos += length
        elif op in ("D", "N"):
            pos += length
    return intervals


# ── Interval / rate helpers ──────────────────────────────────────────────────

def aligned_bp(aligned_intervals, query_intervals):
    total = 0
    i = j = 0
    while i < len(aligned_intervals) and j < len(query_intervals):
        a_start, a_end = aligned_intervals[i]
        q_start, q_end = query_intervals[j]
        overlap = min(a_end, q_end) - max(a_start, q_start)
        if overlap > 0:
            total += overlap
        if a_end < q_end:
            i += 1
        else:
            j += 1
    return total

def count_indel_bases_in_intervals(indel_intervals, query_intervals):
    """
    Count total bases of indel intervals overlapping query intervals.
    Both must be sorted lists of (start, end).
    Uses a two-pointer sweep (same logic as aligned_bp).
    """
    total = 0
    i = j = 0
    while i < len(indel_intervals) and j < len(query_intervals):
        a_start, a_end = indel_intervals[i]
        q_start, q_end = query_intervals[j]
        overlap = min(a_end, q_end) - max(a_start, q_start)
        if overlap > 0:
            total += overlap
        if a_end < q_end:
            i += 1
        else:
            j += 1
    return total

def intervals_overlap(a, b):
    return not (a[1] <= b[0] or a[0] >= b[1])

def random_intervals_like(real_intervals, bounds, rng, max_tries=1000, sample_info=None):
    fake_intervals = []
    lengths = [end - start for start, end in real_intervals]
    for L in lengths:
        start_min = bounds[0]
        start_max = bounds[1] - L
        if start_max < start_min:
            msg = f"Interval length {L} exceeds ASAT bounds {bounds}"
            if sample_info:
                msg += f" for {sample_info}, real intervals: {real_intervals}"
            raise ValueError(msg)
        for _ in range(max_tries):
            s = rng.integers(start_min, start_max + 1)
            candidate = (s, s + L)
            if any(intervals_overlap(candidate, x) for x in fake_intervals):
                continue
            fake_intervals.append(candidate)
            break
        else:
            msg = f"Failed to place non-overlapping fake intervals after {max_tries} tries"
            if sample_info:
                msg += f" for {sample_info}, real intervals: {real_intervals}, bounds: {bounds}"
            raise RuntimeError(msg)
    return sorted(fake_intervals)

def variant_rate(indel_intervals, aligned_intervals, query_intervals):
    aligned = aligned_bp(aligned_intervals, query_intervals)
    if aligned == 0:
        return np.nan
    n = count_indel_bases_in_intervals(indel_intervals, query_intervals)
    return n / aligned

def complement_intervals(intervals, bounds):
    outside = []
    current = bounds[0]
    for start, end in sorted(intervals):
        if start > current:
            outside.append((current, start))
        current = max(current, end)
    if current < bounds[1]:
        outside.append((current, bounds[1]))
    return outside

def interval_union(intervals_a, intervals_b):
    """Return the merged union of two sets of intervals."""
    all_intervals = sorted(intervals_a + intervals_b)
    if not all_intervals:
        return []
    merged = [all_intervals[0]]
    for start, end in all_intervals[1:]:
        if start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged

def interval_subtract(intervals_a, intervals_b):
    """Return parts of intervals_a not overlapping intervals_b (A - B)."""
    result = []
    for a_start, a_end in sorted(intervals_a):
        pieces = [(a_start, a_end)]
        for b_start, b_end in sorted(intervals_b):
            new_pieces = []
            for p_start, p_end in pieces:
                if b_end <= p_start or b_start >= p_end:
                    new_pieces.append((p_start, p_end))
                else:
                    if p_start < b_start:
                        new_pieces.append((p_start, b_start))
                    if p_end > b_end:
                        new_pieces.append((b_end, p_end))
            pieces = new_pieces
        result.extend(pieces)
    return sorted(result)


# ── Data loading ─────────────────────────────────────────────────────────────

def load_all_data(pairwise_dist_dir, indel_parquet, asat_bed_dir, cdr_idx_csv,
                  local_id_bed_dir, cigar_csv_dir):
    """Load all input data and build lookup dictionaries."""
    print("Loading pairwise distances...")
    files = glob.glob(os.path.join(pairwise_dist_dir, "*_r2_QC_v2_centrolign_pairwise_distance.csv"))
    all_dfs = []
    for f in files:
        basename = os.path.basename(f)
        chr_val = basename.split("_")[0]
        df = pd.read_csv(f, header=None, names=["sample1", "sample2", "direct_pairwise_dist", "chr"])
        df["chr"] = chr_val
        all_dfs.append(df)
    all_pairs_dist_df = pd.concat(all_dfs, ignore_index=True)
    all_pairs_dist_df = all_pairs_dist_df[all_pairs_dist_df["direct_pairwise_dist"] < 0.2]
    all_pairs_dist_df["sample_pair"] = (
        all_pairs_dist_df[["sample1", "sample2"]]
        .apply(lambda x: "_".join(sorted(x)), axis=1)
    )
    pairs_df = all_pairs_dist_df[["sample1", "sample2", "sample_pair", "chr"]].copy()

    print("Loading short indels...")
    indels_raw_df = pd.read_parquet(indel_parquet)

    print("Loading ASAT beds...")
    ASAT_DIR = Path(asat_bed_dir)
    dfs = []
    for bed_file in ASAT_DIR.glob("*_asat_arrays.bed"):
        sample = bed_file.name.replace("_asat_arrays.bed", "")
        df = pd.read_csv(bed_file, sep="\t", header=None, names=["contig", "start", "end", "chr"])
        df["sample"] = sample
        dfs.append(df)
    asat_df = pd.concat(dfs, ignore_index=True)

    # ── CDR intervals from centrodip beds ──
    print("Loading CDR beds from index...")
    cdr_idx_df = pd.read_csv(cdr_idx_csv, usecols=["sample_id", "haplotype", "centrodip_final"])
    cdr_idx_df["sample"] = cdr_idx_df["sample_id"] + "." + cdr_idx_df["haplotype"].astype(str)

    contig_to_sample_chr = {}
    for row in asat_df.itertuples(index=False):
        contig_to_sample_chr[row.contig] = (row.sample, row.chr)

    cdr_dfs = []
    missing = 0
    for row in cdr_idx_df.itertuples(index=False):
        bed_path = row.centrodip_final
        sample = row.sample
        if not os.path.exists(bed_path):
            missing += 1
            continue
        bed = pd.read_csv(bed_path, sep="\t", header=None, usecols=[0, 1, 2],
                           names=["contig", "start", "end"])
        bed = bed[bed["contig"].isin(contig_to_sample_chr)]
        bed["sample"] = bed["contig"].map(lambda c: contig_to_sample_chr[c][0])
        bed["chr"] = bed["contig"].map(lambda c: contig_to_sample_chr[c][1])
        bed = bed[bed["sample"] == sample]
        if len(bed) > 0:
            cdr_dfs.append(bed)
    cdr_df = pd.concat(cdr_dfs, ignore_index=True)
    print(f"  Missing bed files: {missing}")
    print(f"  CDR entries: {len(cdr_df)}, samples: {cdr_df['sample'].nunique()}")

    # ── Top 10% local identity intervals ──
    print("Loading local identity beds...")
    base_dir = Path(local_id_bed_dir)
    dfs = []
    for chr_dir in base_dir.iterdir():
        if not chr_dir.is_dir():
            continue
        chr_name = chr_dir.name
        for bed_file in chr_dir.glob("*_local_id_asat_coords.bed"):
            sample = bed_file.name.split("_")[0]
            df = pd.read_csv(bed_file, sep="\t", header=None, usecols=[0, 1, 2, 4],
                             names=["contig", "start", "end", "local_id"])
            df["sample"] = sample
            df["chr"] = chr_name
            dfs.append(df)
    local_id_df = pd.concat(dfs, ignore_index=True)

    thresholds = (
        local_id_df.groupby(["sample", "chr"])["local_id"]
        .quantile(0.9).rename("threshold")
    )
    top10_df = (
        local_id_df.join(thresholds, on=["sample", "chr"])
        .query("local_id >= threshold")
        .drop(columns="threshold")
    )

    # ── Build lookup dicts ──
    print("Building lookup dictionaries...")

    real_cdr_intervals = defaultdict(list)
    for row in cdr_df.itertuples(index=False):
        contig, start, end, sample, chr_ = row
        real_cdr_intervals[(sample, chr_)].append((start, end))
    for key in real_cdr_intervals:
        real_cdr_intervals[key].sort()
    real_cdr_intervals = dict(real_cdr_intervals)

    real_local_id_intervals = defaultdict(list)
    for row in top10_df.itertuples(index=False):
        contig, start, end, local_id, sample, chr_ = row
        real_local_id_intervals[(sample, chr_)].append((start, end))
    for key in real_local_id_intervals:
        real_local_id_intervals[key].sort()
    real_local_id_intervals = dict(real_local_id_intervals)

    asat_bounds = {}
    for row in asat_df.itertuples(index=False):
        contig, start, end, chr_, sample = row
        asat_bounds[(sample, chr_)] = (start, end)

    # Filter out pairs where either sample is missing CDR or local ID
    print("Filtering pairs missing CDR or local ID entries...")
    def pair_missing_data(row):
        s1, s2, pair, chr_ = row
        for s in (s1, s2):
            if (
                (s, chr_) not in real_cdr_intervals
                or len(real_cdr_intervals[(s, chr_)]) == 0
                or (s, chr_) not in real_local_id_intervals
                or len(real_local_id_intervals[(s, chr_)]) == 0
            ):
                return True
        return False

    n_before = len(pairs_df)
    pairs_df = pairs_df[~pairs_df.apply(pair_missing_data, axis=1)].copy()
    print(f"  Pairs before: {n_before}, after filtering: {len(pairs_df)}")

    # Filter out samples where CDRs > 20% of array size
    print("Filtering dense CDR samples (>20% of array)...")
    dense_samples_chrs = set()
    for key in real_cdr_intervals:
        if key not in asat_bounds:
            continue
        sample, chr_ = key
        cdrs = real_cdr_intervals[key]
        bounds_start, bounds_end = asat_bounds[key]
        total_cdr_bases = sum(end - start for start, end in cdrs)
        total_bounds_bases = bounds_end - bounds_start
        if total_cdr_bases / total_bounds_bases > 0.2:
            dense_samples_chrs.add((sample, chr_))
    print(f"  Number of dense sample/chrs: {len(dense_samples_chrs)}")

    def pair_is_dense(row):
        s1, s2, pair, chr_ = row
        return (s1, chr_) in dense_samples_chrs or (s2, chr_) in dense_samples_chrs

    n_before = len(pairs_df)
    pairs_df = pairs_df[~pairs_df.apply(pair_is_dense, axis=1)].copy()
    print(f"  Original pairs: {n_before}, filtered pairs: {len(pairs_df)}")

    indel_pos = defaultdict(list)
    for row in indels_raw_df.itertuples(index=False):
        indel_pos[(row.sample1, row.chr, row.sample_pair)].append((row.ref_start, row.ref_end))
        indel_pos[(row.sample2, row.chr, row.sample_pair)].append((row.query_start, row.query_end))
    indel_pos = {k: sorted(v) for k, v in indel_pos.items()}

    # Contig start map for CIGAR parsing
    contig_start_map = {}
    for row in asat_df.itertuples(index=False):
        contig_start_map[(row.sample, row.contig, row.chr)] = row.start

    # Aligned intervals from CIGAR strings
    print("Parsing CIGAR strings...")
    aligned_intervals = {}
    chroms = [f.replace(".contig_maps.csv", "")
              for f in os.listdir(cigar_csv_dir) if f.endswith(".contig_maps.csv")]
    for chr_ in chroms:
        csv_path = os.path.join(cigar_csv_dir, f"{chr_}.contig_maps.csv")
        df = pd.read_csv(csv_path)
        for row in df.itertuples(index=False):
            s1, c1, s2, c2, cigar_path = row
            if not os.path.exists(cigar_path):
                continue
            with open(cigar_path) as fh:
                cigar_string = fh.readline().strip()
            parsed = parse_cigar(cigar_string)
            ref_start = contig_start_map.get((s1, c1, chr_))
            query_start = contig_start_map.get((s2, c2, chr_))
            if ref_start is None or query_start is None:
                continue
            ref_intervals = build_aligned_intervals_from_ops(parsed, ref_start)
            query_intervals = build_aligned_intervals_from_ops(reverse_cigar(parsed), query_start)
            pair_key = "_".join(sorted([s1, s2]))
            aligned_intervals[(s1, chr_, pair_key)] = ref_intervals
            aligned_intervals[(s2, chr_, pair_key)] = query_intervals

    print("Data loading complete.")
    return (pairs_df, real_cdr_intervals, real_local_id_intervals,
            asat_bounds, indel_pos, aligned_intervals)


# ── Per-chromosome worker ────────────────────────────────────────────────────

_pairs_df = None
_real_cdr_intervals = None
_real_local_id_intervals = None
_asat_bounds = None
_indel_pos = None
_aligned_intervals = None
_n_reps = None

def _init_worker(pairs_df, real_cdr_intervals, real_local_id_intervals,
                 asat_bounds, indel_pos, aligned_intervals, n_reps):
    global _pairs_df, _real_cdr_intervals, _real_local_id_intervals
    global _asat_bounds, _indel_pos, _aligned_intervals, _n_reps
    _pairs_df = pairs_df
    _real_cdr_intervals = real_cdr_intervals
    _real_local_id_intervals = real_local_id_intervals
    _asat_bounds = asat_bounds
    _indel_pos = indel_pos
    _aligned_intervals = aligned_intervals
    _n_reps = n_reps

def _run_chromosome(chr_):
    """Run the permutation test for a single chromosome."""
    print(f"  Starting {chr_} ({_n_reps} reps)...")

    chr_pairs = _pairs_df[_pairs_df['chr'] == chr_]

    pair_rngs = {}
    for idx, r in enumerate(chr_pairs.itertuples(index=False)):
        s1, s2, pair, _ = r
        pair_rngs[pair] = np.random.default_rng(42 + idx)

    null_distr = np.empty(_n_reps)

    for rep in range(_n_reps):
        fake_in_rates = []
        fake_out_rates = []

        for r in chr_pairs.itertuples(index=False):
            s1, s2, pair, _ = r

            for sample in (s1, s2):
                key = (sample, chr_)
                full_key = (sample, chr_, pair)

                if (
                    key not in _asat_bounds
                    or key not in _real_cdr_intervals
                    or key not in _real_local_id_intervals
                    or full_key not in _indel_pos
                    or full_key not in _aligned_intervals
                ):
                    continue

                asat_start, asat_end = _asat_bounds[key]
                rng = pair_rngs[pair]

                # Place fake CDR and fake top 10% independently
                fake_cdr = random_intervals_like(
                    real_intervals=_real_cdr_intervals[key],
                    bounds=(asat_start, asat_end),
                    rng=rng,
                    sample_info=f"{sample}, {chr_}, {pair}, CDR"
                )
                fake_top10 = random_intervals_like(
                    real_intervals=_real_local_id_intervals[key],
                    bounds=(asat_start, asat_end),
                    rng=rng,
                    sample_info=f"{sample}, {chr_}, {pair}, top10"
                )

                # test_region = fake_top10 - fake_CDR
                test_region = interval_subtract(fake_top10, fake_cdr)
                # comparison_region = complement(fake_CDR ∪ fake_top10)
                union_region = interval_union(fake_cdr, fake_top10)
                comparison_region = complement_intervals(union_region, (asat_start, asat_end))

                pos = _indel_pos[full_key]
                aln = _aligned_intervals[full_key]

                rate_in = variant_rate(pos, aln, test_region)
                rate_out = variant_rate(pos, aln, comparison_region)

                if np.isfinite(rate_in) and np.isfinite(rate_out):
                    fake_in_rates.append(rate_in)
                    fake_out_rates.append(rate_out)

        if fake_in_rates and fake_out_rates:
            null_distr[rep] = np.mean(fake_in_rates) - np.mean(fake_out_rates)
        else:
            null_distr[rep] = np.nan

    print(f"  Finished {chr_}")
    return chr_, null_distr


# ── Observed test statistic ──────────────────────────────────────────────────

def compute_observed_test_stats(pairs_df, real_cdr_intervals, real_local_id_intervals,
                                asat_bounds, indel_pos, aligned_intervals):
    """Compute the observed test statistic for each chromosome."""
    results = []
    for chrom in pairs_df['chr'].unique():
        chr_pairs = pairs_df[pairs_df['chr'] == chrom]
        real_in_rates = []
        real_out_rates = []

        for r in chr_pairs.itertuples(index=False):
            s1, s2, pair, chr_ = r

            for sample in (s1, s2):
                key = (sample, chr_)
                full_key = (sample, chr_, pair)
                if (
                    key not in asat_bounds
                    or key not in real_cdr_intervals
                    or key not in real_local_id_intervals
                    or full_key not in indel_pos
                    or full_key not in aligned_intervals
                ):
                    continue

                cdrs = real_cdr_intervals[key]
                top10 = real_local_id_intervals[key]
                asat_start, asat_end = asat_bounds[key]

                # test_region = real top 10% - real CDR
                test_region = interval_subtract(top10, cdrs)
                # comparison_region = complement(real CDR ∪ real top 10%)
                union_region = interval_union(cdrs, top10)
                comparison_region = complement_intervals(union_region, (asat_start, asat_end))

                pos = indel_pos[full_key]
                aln = aligned_intervals[full_key]

                rate_in = variant_rate(pos, aln, test_region)
                rate_out = variant_rate(pos, aln, comparison_region)

                if np.isfinite(rate_in) and np.isfinite(rate_out):
                    real_in_rates.append(rate_in)
                    real_out_rates.append(rate_out)

        if real_in_rates and real_out_rates:
            test_stat = np.mean(real_in_rates) - np.mean(real_out_rates)
        else:
            test_stat = np.nan

        results.append({
            "chr": chrom,
            "n_pairs": len(chr_pairs),
            "n_rates": len(real_in_rates),
            "mean_in_rate": np.mean(real_in_rates) if real_in_rates else np.nan,
            "mean_out_rate": np.mean(real_out_rates) if real_out_rates else np.nan,
            "test_stat": test_stat,
        })

    return pd.DataFrame(results)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Per-chromosome localID-minus-CDR permutation test for short indels (parallelized)"
    )
    parser.add_argument("--n-reps", type=int, default=1000, help="Number of permutation reps (default: 1000)")
    parser.add_argument("--n-workers", type=int, default=8, help="Number of parallel workers (default: 8)")
    parser.add_argument("--output-dir", type=str, required=True, help="Directory to write output CSVs")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Input paths
    pairwise_dist_dir = "/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/"
    indel_parquet = "/private/groups/patenlab/mira/centrolign/analysis/short_indels_pairwise_asm_coords/short_indels_raw_df_tri_traps_distlt0.2.parquet"
    asat_bed_dir = "/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds"
    cdr_idx_csv = "/private/groups/migalab/jmmenend/HPRC/cenSatProject/CDR_data/hprc_ont_centrodip_output_data_index.csv"
    local_id_bed_dir = "/private/groups/patenlab/mira/centrolign/analysis/local_identity/permutation/local_id_beds_asm_coords"
    cigar_csv_dir = "/private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/aligned_bases_per_bed/dist_0.2_smp_contig_maps_cigars"

    # Load data
    (pairs_df, real_cdr_intervals, real_local_id_intervals,
     asat_bounds, indel_pos, aligned_intervals) = load_all_data(
        pairwise_dist_dir, indel_parquet, asat_bed_dir, cdr_idx_csv,
        local_id_bed_dir, cigar_csv_dir
    )

    chromosomes = sorted(pairs_df['chr'].unique())
    print(f"\nRunning {args.n_reps} permutations across {len(chromosomes)} chromosomes "
          f"with {args.n_workers} workers...\n")

    # Parallel permutation tests
    with Pool(
        processes=args.n_workers,
        initializer=_init_worker,
        initargs=(pairs_df, real_cdr_intervals, real_local_id_intervals,
                  asat_bounds, indel_pos, aligned_intervals, args.n_reps),
    ) as pool:
        results = pool.map(_run_chromosome, chromosomes)

    # Collect null distributions
    nulls_per_chr = {chr_: null_distr for chr_, null_distr in results}
    null_df = pd.DataFrame(nulls_per_chr)
    null_out = os.path.join(args.output_dir, "localid_minus_cdr_indel_null_distribution_per_chr.csv")
    null_df.to_csv(null_out, index=False)
    print(f"\nNull distributions saved to {null_out}")

    # Compute observed test statistics
    print("\nComputing observed test statistics...")
    per_chrom_results_df = compute_observed_test_stats(
        pairs_df, real_cdr_intervals, real_local_id_intervals,
        asat_bounds, indel_pos, aligned_intervals
    )

    # Compute two-sided p-values
    per_chrom_results_df["p_value"] = np.nan
    for _, row in per_chrom_results_df.iterrows():
        chrom = row["chr"]
        test_stat = row["test_stat"]
        if chrom not in null_df.columns or np.isnan(test_stat):
            continue
        null_vals = null_df[chrom].dropna().values
        p_value = np.mean(np.abs(null_vals) >= abs(test_stat))
        if p_value == 0:
            p_value = 1 / len(null_vals)
        per_chrom_results_df.loc[per_chrom_results_df["chr"] == chrom, "p_value"] = p_value

    results_out = os.path.join(args.output_dir, "localid_minus_cdr_indel_per_chrom_results.csv")
    per_chrom_results_df.to_csv(results_out, index=False)
    print(f"Per-chrom results saved to {results_out}")

    # Print summary
    print("\n" + "=" * 80)
    print("Results summary (top10% local ID minus CDR, short indels):")
    print("=" * 80)
    for _, row in per_chrom_results_df.iterrows():
        p = row['p_value']
        p_str = f"{p:.2e}" if p < 1e-5 else f"{p:.6f}"
        print(f"  {row['chr']:6s}  test_stat: {row['test_stat']:+.6f}  p_value: {p_str}")


if __name__ == "__main__":
    main()
