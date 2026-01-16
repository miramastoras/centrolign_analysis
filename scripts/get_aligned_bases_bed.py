import re
import csv
import os
import argparse

import os
import csv
import re

# -------------------------
# CIGAR helpers
# -------------------------
def parse_cigar(cigar):
    """Parse CIGAR string into [(length, op), ...]"""
    return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]

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
    """
    Build aligned intervals (start, end) along coordinate system starting at start_pos.
    """
    intervals = []
    pos = start_pos
    for length, op in cigar_ops:
        if op in ("M", "=", "X"):
            intervals.append((pos, pos + length))
            pos += length
        elif op in ("D", "N"):  # Consume reference only
            pos += length
        else:  # I, S, H, P
            # I consumes query in query coordinates
            if op == "I":
                pos += length
            # S, H, P do not consume anything
    return intervals

def count_aligned_bases(intervals, start, end):
    count = 0
    for a_start, a_end in intervals:
        count += max(0, min(end, a_end) - max(start, a_start))
    return count

# -------------------------
# BED helpers
# -------------------------
def read_bed_subset(bed_file, contig):
    """
    Read BED file, subset by contig.
    Returns:
        windows: list of full BED lines as lists of columns
        min_start: minimum start coordinate for that contig
    """
    windows = []
    min_start = None
    with open(bed_file) as f:
        for line in f:
            if line.strip() == "" or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if chrom != contig:
                continue
            windows.append(fields)  # Keep the full BED row as a list
            if min_start is None or start < min_start:
                min_start = start
    return windows, min_start


# -------------------------
# CIGAR indexing
# -------------------------
def index_cigar_files(cigar_dir):
    """
    Build dictionary: frozenset({sample1, sample2}) -> (ref_sample, query_sample, filepath)
    Keeps the first file encountered for duplicate pairs
    """
    index = {}
    for fname in sorted(os.listdir(cigar_dir)):
        if not fname.startswith("pairwise_cigar_") or not fname.endswith(".txt"):
            continue
        parts = fname.rstrip(".txt").split("_")
        if len(parts) != 4:
            continue
        _, _, ref, query = parts
        key = frozenset((ref, query))
        if key not in index:
            index[key] = (ref, query, os.path.join(cigar_dir, fname))
    return index

def index_bed_files(bed_dir, suffix):
    """
    Build dictionary:
      frozenset({sample1, sample2}) -> filepath
    Keeps the first BED file found for duplicates.
    """
    bed_index = {}
    for fname in sorted(os.listdir(bed_dir)):
        if not fname.endswith(f".{suffix}.bed"):
            continue
        base = fname.rsplit(f".{suffix}.bed", 1)[0]
        parts = base.split("_")
        if len(parts) != 2:
            continue
        s1, s2 = parts
        key = frozenset((s1, s2))
        if key not in bed_index:
            bed_index[key] = os.path.join(bed_dir, fname)
    return bed_index


# -------------------------
# Main wrapper
# -------------------------
def run_cigar_bed_wrapper(csv_file, cigar_dir, bed_dir, bed_suffix, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    cigar_index = index_cigar_files(cigar_dir)
    bed_index = index_bed_files(bed_dir,bed_suffix)

    with open(csv_file) as csvfh:
        reader = csv.reader(csvfh)
        for row in reader:
            # Columns: sample1, sample1_contig, sample2, sample2_contig
            s1 = row[0]
            c1 = row[1]
            s2 = row[2]
            c2 = row[3]

            key = frozenset((s1, s2))
            if key not in cigar_index:
                raise FileNotFoundError(f"No cigar file for {s1}, {s2}")

            ref_sample, query_sample, cigar_file = cigar_index[key]

            # Map contigs to ref/query
            ref_contig = c1 if ref_sample == s1 else c2
            query_contig = c2 if query_sample == s2 else c1

            if key not in bed_index:
                raise FileNotFoundError(f"No BED file found for pair {s1}, {s2}")
            
            bed_file = bed_index[key]

            ref_windows, ref_start = read_bed_subset(bed_file, ref_contig)
            query_windows, query_start = read_bed_subset(bed_file, query_contig)

            if ref_start is None or query_start is None:
                raise ValueError(f"No BED entries for contigs {ref_contig} / {query_contig}")

            # Read cigar
            with open(cigar_file) as f:
                cigar = f.readline().strip()
            parsed_cigar = parse_cigar(cigar)

            # Build aligned intervals
            ref_intervals = build_aligned_intervals_from_ops(parsed_cigar, ref_start)
            query_intervals = build_aligned_intervals_from_ops(
                reverse_cigar(parsed_cigar), query_start
            )

            # Output
            out_path = os.path.join(output_dir, f"{s1}_{s2}.aligned_bases.bed")

            with open(out_path, "w") as out:
                for fields in ref_windows:
                    chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                    ref_count = count_aligned_bases(ref_intervals, start, end)
                    out.write("\t".join(fields) + f"\t{ref_count}\n")
                for fields in query_windows:
                    chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                    query_count = count_aligned_bases(query_intervals, start, end)
                    out.write("\t".join(fields) + f"\t{query_count}\n")

            print(f"Processed {s1} vs {s2}")


import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Report aligned bases per BED window from pairwise CIGAR string."
    )
    parser.add_argument("-c", "--csv", required=True, help="CSV file (no header) with sample pairs")
    parser.add_argument("-b", "--bed_dir", required=True, help="Directory containing BED files. Assumes both ref and query windows are in the same bed file")
    parser.add_argument("-s", "--bed_suffix", required=True, help="Suffix used in BED files")
    parser.add_argument("-g", "--cigar_dir", required=True, help="Directory containing CIGAR files")
    parser.add_argument("-o", "--out_dir", required=True, help="Directory to write output files")

    args = parser.parse_args()

    run_cigar_bed_wrapper(args.csv, args.cigar_dir, args.bed_dir, args.bed_suffix, args.out_dir)


if __name__ == "__main__":
    main()
