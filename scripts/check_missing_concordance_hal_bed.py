#!/usr/bin/env python3
"""
For each sample in a VCF, check that missing alleles (.\/. or .\/X) correspond
to positions NOT covered by that haplotype's BED file.

Haplotype BED files are named {sample}.1.bed and {sample}.2.bed, matching
allele 1 and allele 2 of the diploid genotype respectively.

Per-sample summary output (TSV):
  sample,
  total_variants,
  n_missing_alleles,         -- total missing alleles (each GT has 2 alleles)
  n_missing_and_uncovered,   -- missing allele AND haplotype BED has no coverage (expected)
  n_missing_and_covered,     -- missing allele BUT haplotype BED has coverage (suspicious)
  n_present_and_uncovered,   -- allele called BUT haplotype BED has no coverage (suspicious)
  n_present_and_covered,     -- allele called AND haplotype BED has coverage (expected)
  pct_missing_explained      -- n_missing_and_uncovered / n_missing_alleles * 100

Usage:
    python check_missing_concordance.py \
        --vcf chr13_flanking_50kb.vcf.gz \
        --beddir per_sample_beds/ \
        --out chr13_concordance.tsv
"""

import os
import gzip
import argparse
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--vcf",    required=True, help="VCF or VCF.gz")
    p.add_argument("--beddir", required=True, help="Directory of per-sample haplotype BED files")
    p.add_argument("--out",    required=True, help="Output TSV")
    return p.parse_args()


def load_bed_coverage(bed_path):
    """
    Load a BED file into dict of chrom -> sorted list of (start, end).
    Returns empty dict if file doesn't exist.
    """
    cov = defaultdict(list)
    if not os.path.exists(bed_path):
        return cov
    with open(bed_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            cov[chrom].append((start, end))
    for chrom in cov:
        cov[chrom].sort()
    return cov


def is_covered(cov, chrom, pos0):
    """Binary search: is pos (0-based) covered by any interval in cov[chrom]?"""
    intervals = cov.get(chrom, [])
    lo, hi = 0, len(intervals) - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        start, end = intervals[mid]
        if pos0 < start:
            hi = mid - 1
        elif pos0 >= end:
            lo = mid + 1
        else:
            return True
    return False


def open_vcf(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def main():
    args = parse_args()

    # --- Read VCF header to get sample names ---
    samples = []
    with open_vcf(args.vcf) as f:
        for line in f:
            if line.startswith("#CHROM"):
                samples = line.strip().split("\t")[9:]
                break

    print(f"Found {len(samples)} samples in VCF")

    # --- Load haplotype BED coverage per sample ---
    print("Loading BED files...")
    # cov[sample][hap] where hap is 0 (allele 1) or 1 (allele 2)
    cov = {}
    missing_beds = []
    for samp in samples:
        hap1 = load_bed_coverage(os.path.join(args.beddir, f"{samp}.1.bed"))
        hap2 = load_bed_coverage(os.path.join(args.beddir, f"{samp}.2.bed"))
        cov[samp] = [hap1, hap2]
        if not hap1 and not hap2:
            missing_beds.append(samp)

    if missing_beds:
        print(f"  WARNING: {len(missing_beds)} samples have no BED files at all")

    # --- Per-sample counters ---
    # [total_variants, n_missing_alleles, missing_uncov, missing_cov, present_uncov, present_cov]
    counts = {s: [0, 0, 0, 0, 0, 0] for s in samples}

    # --- Stream VCF ---
    print("Streaming VCF...")
    with open_vcf(args.vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            pos0  = int(fields[1]) - 1  # 1-based -> 0-based
            fmt   = fields[8].split(":")
            gt_idx = fmt.index("GT") if "GT" in fmt else 0

            for i, samp in enumerate(samples):
                sample_field = fields[9 + i].split(":")
                gt = sample_field[gt_idx] if gt_idx < len(sample_field) else "./."
                alleles = gt.replace("|", "/").split("/")

                c = counts[samp]
                c[0] += 1  # total_variants

                for hap_idx, allele in enumerate(alleles):
                    missing = allele == "."
                    covered = is_covered(cov[samp][hap_idx], chrom, pos0)

                    if missing:
                        c[1] += 1  # n_missing_alleles
                        if not covered:
                            c[2] += 1  # missing_and_uncovered (expected)
                        else:
                            c[3] += 1  # missing_and_covered   (suspicious)
                    else:
                        if not covered:
                            c[4] += 1  # present_and_uncovered (suspicious)
                        else:
                            c[5] += 1  # present_and_covered   (expected)

    # --- Write output ---
    print(f"Writing {args.out}...")
    with open(args.out, "w") as out:
        out.write(
            "sample\ttotal_variants\tn_missing_alleles\tn_missing_and_uncovered\t"
            "n_missing_and_covered\tn_present_and_uncovered\tn_present_and_covered\t"
            "pct_missing_explained\n"
        )
        for samp in samples:
            c = counts[samp]
            total, n_miss, miss_uncov, miss_cov, pres_uncov, pres_cov = c
            pct = 100 * miss_uncov / n_miss if n_miss > 0 else float("nan")
            out.write(
                f"{samp}\t{total}\t{n_miss}\t{miss_uncov}\t"
                f"{miss_cov}\t{pres_uncov}\t{pres_cov}\t{pct:.1f}\n"
            )

    print("Done.")


if __name__ == "__main__":
    main()