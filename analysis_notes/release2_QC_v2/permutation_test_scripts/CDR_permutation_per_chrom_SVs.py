#!/usr/bin/env python3
"""
Per-chromosome permutation test: do CDRs
have an elevated SV rate compared to the rest of the array?

Parallelizes across chromosomes using multiprocessing.

Usage:
    python CDR_permutation_per_chrom_SVs.py --n-reps 1000 --n-workers 8 --output-dir /path/to/output
"""

import argparse
import os
import glob
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool

import numpy as np
import pandas as pd


# ── Interval / rate helpers ──────────────────────────────────────────────────

def count_sv_bases_in_intervals(sv_intervals, query_intervals):
    """
    Count total bases of SV intervals overlapping query intervals.
    Both must be sorted lists of (start, end).
    Uses a two-pointer sweep.
    """
    total = 0
    i = j = 0
    while i < len(sv_intervals) and j < len(query_intervals):
        a_start, a_end = sv_intervals[i]
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
            msg = f"CDR length {L} exceeds ASAT bounds {bounds}"
            if sample_info:
                msg += f" for {sample_info}, real CDRs: {real_intervals}"
            raise ValueError(msg)
        for _ in range(max_tries):
            s = rng.integers(start_min, start_max + 1)
            candidate = (s, s + L)
            if any(intervals_overlap(candidate, x) for x in fake_intervals):
                continue
            fake_intervals.append(candidate)
            break
        else:
            msg = f"Failed to place non-overlapping fake CDRs after {max_tries} tries"
            if sample_info:
                msg += f" for {sample_info}, real CDRs: {real_intervals}, bounds: {bounds}"
            raise RuntimeError(msg)
    return sorted(fake_intervals)

def variant_rate(sv_intervals, query_intervals):
    window_size = sum(end - start for start, end in query_intervals)
    if window_size == 0:
        return np.nan
    n = count_sv_bases_in_intervals(sv_intervals, query_intervals)
    return n / window_size

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


# ── Data loading ─────────────────────────────────────────────────────────────

def load_all_data(pairwise_dist_dir, sv_parquet, asat_bed_dir, cdr_idx_csv):
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

    print("Loading SVs...")
    svs_raw_df = pd.read_parquet(sv_parquet)

    print("Loading ASAT beds...")
    ASAT_DIR = Path(asat_bed_dir)
    dfs = []
    for bed_file in ASAT_DIR.glob("*_asat_arrays.bed"):
        sample = bed_file.name.replace("_asat_arrays.bed", "")
        df = pd.read_csv(bed_file, sep="\t", header=None, names=["contig", "start", "end", "chr"])
        df["sample"] = sample
        dfs.append(df)
    asat_df = pd.concat(dfs, ignore_index=True)

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

    # ── Build lookup dicts ──
    print("Building lookup dictionaries...")

    real_cdr_intervals = defaultdict(list)
    for row in cdr_df.itertuples(index=False):
        contig, start, end, sample, chr_ = row
        real_cdr_intervals[(sample, chr_)].append((start, end))
    for key in real_cdr_intervals:
        real_cdr_intervals[key].sort()
    real_cdr_intervals = dict(real_cdr_intervals)

    asat_bounds = {}
    for row in asat_df.itertuples(index=False):
        contig, start, end, chr_, sample = row
        asat_bounds[(sample, chr_)] = (start, end)

    # Filter out pairs where either sample has no CDR entries
    print("Filtering pairs missing CDR entries...")
    def pair_missing_cdr(row):
        s1, s2, pair, chr_ = row
        return (
            (s1, chr_) not in real_cdr_intervals
            or (s2, chr_) not in real_cdr_intervals
            or len(real_cdr_intervals.get((s1, chr_), [])) == 0
            or len(real_cdr_intervals.get((s2, chr_), [])) == 0
        )
    n_before_cdr = len(pairs_df)
    pairs_df = pairs_df[~pairs_df.apply(pair_missing_cdr, axis=1)].copy()
    print(f"  Pairs before: {n_before_cdr}, after filtering missing CDRs: {len(pairs_df)}")

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

    sv_pos = defaultdict(list)
    for row in svs_raw_df.itertuples(index=False):
        sv_pos[(row.sample1, row.chr, row.sample_pair)].append((row.ref_start, row.ref_end))
        sv_pos[(row.sample2, row.chr, row.sample_pair)].append((row.query_start, row.query_end))
    sv_pos = {k: sorted(v) for k, v in sv_pos.items()}

    print("Data loading complete.")
    return pairs_df, real_cdr_intervals, asat_bounds, sv_pos


# ── Per-chromosome worker ────────────────────────────────────────────────────

# These globals are set once by each worker via the initializer
_pairs_df = None
_real_cdr_intervals = None
_asat_bounds = None
_sv_pos = None
_n_reps = None

def _init_worker(pairs_df, real_cdr_intervals, asat_bounds, sv_pos, n_reps):
    global _pairs_df, _real_cdr_intervals, _asat_bounds, _sv_pos, _n_reps
    _pairs_df = pairs_df
    _real_cdr_intervals = real_cdr_intervals
    _asat_bounds = asat_bounds
    _sv_pos = sv_pos
    _n_reps = n_reps

def _run_chromosome(chr_):
    """Run the permutation test for a single chromosome. Returns (chr_, null_distr)."""
    print(f"  Starting {chr_} ({_n_reps} reps)...")

    chr_pairs = _pairs_df[_pairs_df['chr'] == chr_]

    # Prepare one independent RNG per pair
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

            if (
                (s1, chr_) not in _real_cdr_intervals
                or (s2, chr_) not in _real_cdr_intervals
                or len(_real_cdr_intervals[(s1, chr_)]) == 0
                or len(_real_cdr_intervals[(s2, chr_)]) == 0
            ):
                continue

            for sample in (s1, s2):
                key = (sample, chr_)
                full_key = (sample, chr_, pair)

                if (
                    key not in _asat_bounds
                    or full_key not in _sv_pos
                ):
                    continue

                asat_start, asat_end = _asat_bounds[key]
                real_cdrs = _real_cdr_intervals[key]

                fake_cdr = random_intervals_like(
                    real_intervals=real_cdrs,
                    bounds=(asat_start, asat_end),
                    rng=pair_rngs[pair],
                    sample_info=f"{sample}, {chr_}, {pair}"
                )

                fake_outside = complement_intervals(fake_cdr, (asat_start, asat_end))

                pos = _sv_pos[full_key]

                rate_in = variant_rate(pos, fake_cdr)
                rate_out = variant_rate(pos, fake_outside)

                if np.isfinite(rate_in) and np.isfinite(rate_out):
                    fake_in_rates.append(rate_in)
                    fake_out_rates.append(rate_out)

        null_distr[rep] = np.mean(fake_in_rates) - np.mean(fake_out_rates)

    print(f"  Finished {chr_}")
    return chr_, null_distr


# ── Per-chromosome test statistic (observed) ─────────────────────────────────

def compute_observed_test_stats(pairs_df, real_cdr_intervals, asat_bounds, sv_pos):
    """Compute the observed test statistic for each chromosome."""
    results = []
    for chrom in pairs_df['chr'].unique():
        chr_pairs = pairs_df[pairs_df['chr'] == chrom]
        real_in_rates = []
        real_out_rates = []

        for r in chr_pairs.itertuples(index=False):
            s1, s2, pair, chr_ = r
            if (
                (s1, chr_) not in real_cdr_intervals
                or (s2, chr_) not in real_cdr_intervals
                or len(real_cdr_intervals[(s1, chr_)]) == 0
                or len(real_cdr_intervals[(s2, chr_)]) == 0
            ):
                continue
            for sample in (s1, s2):
                key = (sample, chr_)
                full_key = (sample, chr_, pair)
                if (
                    key not in asat_bounds
                    or full_key not in sv_pos
                ):
                    continue
                cdrs = real_cdr_intervals[key]
                asat_start, asat_end = asat_bounds[key]
                outside = complement_intervals(cdrs, (asat_start, asat_end))
                pos = sv_pos[full_key]
                rate_in = variant_rate(pos, cdrs)
                rate_out = variant_rate(pos, outside)
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
        description="Per-chromosome CDR permutation test for SVs (parallelized)"
    )
    parser.add_argument("--n-reps", type=int, default=1000, help="Number of permutation reps (default: 1000)")
    parser.add_argument("--n-workers", type=int, default=8, help="Number of parallel workers (default: 8)")
    parser.add_argument("--output-dir", type=str, required=True, help="Directory to write output CSVs")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Input paths
    pairwise_dist_dir = "/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/"
    sv_parquet = "/private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise_asm_coords/svs_df_tri_traps_distlt0.2.parquet"
    asat_bed_dir = "/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds"
    cdr_idx_csv = "/private/groups/migalab/jmmenend/HPRC/cenSatProject/CDR_data/hprc_ont_centrodip_output_data_index.csv"

    # Load data (serial, one-time)
    pairs_df, real_cdr_intervals, asat_bounds, sv_pos = load_all_data(
        pairwise_dist_dir, sv_parquet, asat_bed_dir, cdr_idx_csv
    )

    chromosomes = sorted(pairs_df['chr'].unique())
    print(f"\nRunning {args.n_reps} permutations across {len(chromosomes)} chromosomes with {args.n_workers} workers...\n")

    # Parallel permutation tests
    with Pool(
        processes=args.n_workers,
        initializer=_init_worker,
        initargs=(pairs_df, real_cdr_intervals, asat_bounds, sv_pos, args.n_reps),
    ) as pool:
        results = pool.map(_run_chromosome, chromosomes)

    # Collect null distributions
    nulls_per_chr = {chr_: null_distr for chr_, null_distr in results}
    null_df = pd.DataFrame(nulls_per_chr)
    null_out = os.path.join(args.output_dir, "cdr_sv_null_distribution_per_chr.csv")
    null_df.to_csv(null_out, index=False)
    print(f"\nNull distributions saved to {null_out}")

    # Compute observed test statistics
    print("\nComputing observed test statistics...")
    per_chrom_results_df = compute_observed_test_stats(
        pairs_df, real_cdr_intervals, asat_bounds, sv_pos
    )

    # Compute p-values
    per_chrom_results_df["p_value"] = np.nan
    for _, row in per_chrom_results_df.iterrows():
        chrom = row["chr"]
        test_stat = row["test_stat"]
        if chrom not in null_df.columns:
            continue
        null_vals = null_df[chrom].values
        p_value = np.mean(np.abs(null_vals) >= abs(test_stat))
        if p_value == 0:
            p_value = 1 / len(null_vals)
        per_chrom_results_df.loc[per_chrom_results_df["chr"] == chrom, "p_value"] = p_value

    results_out = os.path.join(args.output_dir, "cdr_sv_per_chrom_results.csv")
    per_chrom_results_df.to_csv(results_out, index=False)
    print(f"Per-chrom results saved to {results_out}")

    # Print summary
    print("\n" + "=" * 80)
    print("Results summary (SVs):")
    print("=" * 80)
    for _, row in per_chrom_results_df.iterrows():
        p = row['p_value']
        p_str = f"{p:.2e}" if p < 1e-5 else f"{p:.6f}"
        print(f"  {row['chr']:6s}  test_stat: {row['test_stat']:+.6f}  p_value: {p_str}")


if __name__ == "__main__":
    main()
