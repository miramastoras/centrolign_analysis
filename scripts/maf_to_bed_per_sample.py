#!/usr/bin/env python3
"""
Parse a MAF file and produce one BED file per sample with aligned segments
in CHM13 reference coordinates.

Output BED columns:
  chrom, start, end, sample, aligned_bases, gap_count

Usage:
    python maf_to_bed_per_sample.py \
        --maf input.maf \
        --ref CHM13 \
        --outdir per_sample_beds/
"""

import os
import argparse


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--maf", required=True)
    p.add_argument("--ref", default="CHM13", help="Reference genome name in MAF (default: CHM13)")
    p.add_argument("--outdir", default="per_sample_beds")
    return p.parse_args()


def sample_name(src):
    """Strip genome prefix: 'HG00097.1.CM094066.1' -> 'HG00097.1'"""
    parts = src.split(".")
    return ".".join(parts[:2]) if len(parts) > 1 else src


def chrom_name(src):
    """Strip genome prefix: 'CHM13.chr1' -> 'chr1'"""
    parts = src.split(".")
    return parts[1] if len(parts) > 1 else src


def process_block(block, ref, write_record):
    ref_prefix = ref + "."
    ref_comp = None
    sample_comps = []
    for fields in block:
        src, start, size, strand, src_size, text = (
            fields[1], int(fields[2]), int(fields[3]),
            fields[4], int(fields[5]), fields[6]
        )
        if src.startswith(ref_prefix) or src == ref:
            ref_comp = (src, start, size, strand, src_size, text)
        else:
            sample_comps.append((src, start, size, strand, src_size, text))

    if ref_comp is None:
        return

    ref_src, ref_start, ref_size, _, _, ref_text = ref_comp
    chrom = chrom_name(ref_src)
    ref_bed_start = ref_start
    ref_bed_end   = ref_start + ref_size

    for src, start, size, strand, src_size, text in sample_comps:
        samp = sample_name(src)
        if samp.startswith("_"):
            continue
        gap_count = sum(
            1 for ref_b, b in zip(ref_text, text) if ref_b != "-" and b == "-"
        )
        write_record(samp, f"{chrom}\t{ref_bed_start}\t{ref_bed_end}\t{samp}\t{gap_count}\n")


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    ref_prefix = args.ref + "."
    seen_samples = set()  # track which files have been created (to use "w" vs "a")

    def write_record(samp, record):
        mode = "a" if samp in seen_samples else "w"
        seen_samples.add(samp)
        with open(os.path.join(args.outdir, f"{samp}.bed"), mode) as fh:
            fh.write(record)

    def process_block(block):
        ref_comp = None
        sample_comps = []
        for fields in block:
            src, start, size, strand, src_size, text = (
                fields[1], int(fields[2]), int(fields[3]),
                fields[4], int(fields[5]), fields[6]
            )
            if src.startswith(ref_prefix) or src == args.ref:
                ref_comp = (src, start, size, strand, src_size, text)
            else:
                sample_comps.append((src, start, size, strand, src_size, text))

        if ref_comp is None:
            return

        ref_src, ref_start, ref_size, _, _, ref_text = ref_comp
        chrom = chrom_name(ref_src)
        ref_bed_start = ref_start
        ref_bed_end   = ref_start + ref_size

        for src, start, size, strand, src_size, text in sample_comps:
            samp = sample_name(src)
            gap_count = sum(
                1 for ref_b, b in zip(ref_text, text) if ref_b != "-" and b == "-"
            )
            write_record(samp, f"{chrom}\t{ref_bed_start}\t{ref_bed_end}\t{samp}\t{gap_count}\n")

    block = []
    with open(args.maf) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("a"):
                block = []
            elif line.startswith("s"):
                block.append(line.split())
            elif line.strip() == "" and block:
                process_block(block, args.ref, write_record)
                block = []

    # Handle last block if file doesn't end with a blank line
    if block:
        process_block(block, args.ref, write_record)

    print(f"Done. Wrote {len(seen_samples)} sample BED files to {args.outdir}/")
    print(f"Samples: {sorted(seen_samples)[:5]} ...")


if __name__ == "__main__":
    main()