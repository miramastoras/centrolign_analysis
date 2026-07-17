#!/usr/bin/env python3
"""
Annotates genes with two columns appended to the input TSV:
  between_satellites_100kb : TRUE if any satellite block > 100 kb exists on BOTH
                             sides of the gene, searched only within the censat
                             array bounds from the regions file
  on_acrocentric_short_arm : TRUE if on an acrocentric chrom AND gene lies left of the
                             active_hor array start on that contig

Usage: annotate_satellite_neighbors.py <tsv> <bed> <out> [--regions-dir <dir>]
"""

import sys, os, re

ACROCENTRIC = {'chr13', 'chr14', 'chr15', 'chr21', 'chr22'}

def strip_contig(name):
    """'sample#hap#CM094067.1' -> 'CM094067.1'; plain names pass through."""
    parts = name.split('#')
    return parts[-1] if len(parts) >= 2 else name

def is_active_hor(name):
    return bool(re.match(r'active_hor', name, re.I))

def is_ct(name):
    return bool(re.match(r'ct', name, re.I))

def parse_bed(bed_path):
    """
    Returns:
      blocks    : contig -> sorted list of (start, end, size)
      hor_starts: contig -> min start of any active_hor block
    """
    blocks = {}
    hor_starts = {}
    with open(bed_path) as f:
        for line in f:
            if line.startswith('track'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 4:
                continue
            contig = strip_contig(cols[0])
            s, e   = int(cols[1]), int(cols[2])
            name   = cols[3]
            blocks.setdefault(contig, []).append((s, e, e - s, name))
            if is_active_hor(name):
                if contig not in hor_starts or s < hor_starts[contig]:
                    hor_starts[contig] = s
    for c in blocks:
        blocks[c].sort()
    return blocks, hor_starts

def load_array_bounds(regions_dir, chrom, sample_id, haplotype):
    """
    Look up the censat array region_start / region_end for this sample+hap+chrom
    from the per-chrom regions TSV files.  Returns (array_start, array_end) or
    (None, None) if the file or row is not found.

    Haplotype matching: the regions file stores 1/2; the script also accepts
    'mat' (→2) and 'pat' (→1) from filenames.
    """
    hap_map = {'mat': '2', 'pat': '1'}
    hap_str = hap_map.get(haplotype, haplotype)

    regions_file = os.path.join(regions_dir, f'censat_regions_pass_qc_{chrom}.tsv')
    if not os.path.exists(regions_file):
        return None, None

    with open(regions_file) as f:
        header = f.readline().rstrip('\n').split('\t')
        try:
            si_col  = header.index('sample_id')
            hap_col = header.index('haplotype')
            rs_col  = header.index('region_start')
            re_col  = header.index('region_end')
        except ValueError:
            return None, None
        for line in f:
            cols = line.rstrip('\n').split('\t')
            if cols[si_col] == sample_id and str(cols[hap_col]) == hap_str:
                return int(cols[rs_col]), int(cols[re_col])
    return None, None

def merge_nonct_blocks(contig_blocks):
    """
    Merge consecutive touching/overlapping non-CT blocks into single 'satellite'
    blocks. CT blocks are kept as-is. Returns a new sorted block list.
    """
    if not contig_blocks:
        return []
    merged = []
    cur_s = cur_e = None
    for s, e, size, name in contig_blocks:   # already sorted by start
        if is_ct(name):
            if cur_s is not None:            # flush pending satellite merge
                merged.append((cur_s, cur_e, cur_e - cur_s, 'satellite'))
                cur_s = cur_e = None
            merged.append((s, e, size, name))
        else:
            if cur_s is None:
                cur_s, cur_e = s, e
            elif s <= cur_e:                 # touching or overlapping — extend
                cur_e = max(cur_e, e)
            else:                            # gap between non-CT blocks
                merged.append((cur_s, cur_e, cur_e - cur_s, 'satellite'))
                cur_s, cur_e = s, e
    if cur_s is not None:
        merged.append((cur_s, cur_e, cur_e - cur_s, 'satellite'))
    return merged

def has_large_flanking_blocks(contig_blocks, gene_start, gene_end,
                               min_size=100_000,
                               array_start=None, array_end=None):
    """
    Return (left_ok, right_ok).
    left_ok  : True if ANY *satellite* block with end  <= gene_start has size > min_size
    right_ok : True if ANY *satellite* block with start >= gene_end   has size > min_size
    CT blocks are excluded — only non-CT satellite blocks count.
    If array_start/array_end are given, only blocks overlapping that window are considered,
    preventing distant arm satellites from qualifying.
    """
    left_ok = right_ok = False
    for s, e, size, name in contig_blocks:
        if is_ct(name):
            continue
        if array_start is not None and e < array_start:
            continue
        if array_end is not None and s > array_end:
            continue
        if e <= gene_start and size > min_size:
            left_ok = True
        if s >= gene_end and size > min_size:
            right_ok = True
        if left_ok and right_ok:
            break
    return left_ok, right_ok

def main():
    # Parse args: positional <tsv> <bed> <out>, optional --regions-dir <dir>
    args = sys.argv[1:]
    regions_dir = None
    if '--regions-dir' in args:
        idx = args.index('--regions-dir')
        regions_dir = args[idx + 1]
        args = args[:idx] + args[idx + 2:]
    if len(args) != 3:
        print("Usage: annotate_satellite_neighbors.py <tsv> <bed> <out> [--regions-dir <dir>]",
              file=sys.stderr)
        sys.exit(1)

    tsv_path, bed_path, out_path = args

    # Extract sample_id, haplotype, chrom from filename
    # e.g. HG00097_hap1_chr13_intersect.tsv  or  HG00408_mat_chr1_intersect.tsv
    fname = os.path.basename(tsv_path)
    m = re.match(r'(.+?)_(hap\d+|mat|pat)_(chr\w+)_', fname)
    sample_id  = m.group(1)                                if m else None
    hap_raw    = m.group(2)                                if m else None
    haplotype  = re.sub(r'hap', '', hap_raw) if hap_raw and hap_raw.startswith('hap') else hap_raw
    chrom      = m.group(3)                                if m else None
    is_acrocentric = chrom in ACROCENTRIC if chrom else False

    # Load censat array bounds for this sample if regions dir provided
    array_start = array_end = None
    if regions_dir and sample_id and haplotype and chrom:
        array_start, array_end = load_array_bounds(regions_dir, chrom, sample_id, haplotype)
        if array_start is None:
            print(f"Warning: no region bounds found for {sample_id} hap={haplotype} {chrom}; "
                  f"satellite search unclamped", file=sys.stderr)

    blocks, hor_starts = parse_bed(bed_path)
    blocks = {c: merge_nonct_blocks(b) for c, b in blocks.items()}

    HEADER = ('contig\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\t'
              'attributes\tcontig2\toverlap_bp\t'
              'between_satellites_100kb\ton_acrocentric_short_arm\n')

    n_written = 0
    with open(tsv_path) as fin, open(out_path, 'w') as fout:
        fout.write(HEADER)
        for line in fin:
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 5:
                fout.write(line.rstrip('\n') + '\tFALSE\tFALSE\n')
                continue

            contig     = cols[0]
            gene_start = int(cols[3])
            gene_end   = int(cols[4])

            # ── between_satellites_100kb ──────────────────────────────────────
            # Apply the 100kb flanking check to all genes. Only non-CT satellite
            # blocks count as valid neighbors (CT blocks are excluded inside
            # has_large_flanking_blocks).
            cb = blocks.get(contig, [])
            left_ok, right_ok = has_large_flanking_blocks(
                cb, gene_start, gene_end,
                array_start=array_start, array_end=array_end)
            between = left_ok and right_ok

            # ── on_acrocentric_short_arm ──────────────────────────────────────
            if is_acrocentric and contig in hor_starts:
                acro_arm = gene_end <= hor_starts[contig]
            else:
                acro_arm = False

            fout.write(line.rstrip('\n')
                       + '\t' + ('TRUE' if between  else 'FALSE')
                       + '\t' + ('TRUE' if acro_arm else 'FALSE')
                       + '\n')
            n_written += 1

    print(f"Done: {n_written} genes written to {out_path}")

if __name__ == '__main__':
    main()
