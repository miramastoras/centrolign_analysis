import re
import csv
import os
import argparse

# -------------------------
# CIGAR helpers
# -------------------------
def parse_cigar(cigar):
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
    intervals = []
    pos = start_pos
    for length, op in cigar_ops:
        if op in ("M", "=", "X"):
            intervals.append((pos, pos + length))
            pos += length
        elif op in ("D", "N"):
            pos += length
        #elif op == "I":
            #pos += length
    return intervals

def count_aligned_bases(intervals, start, end):
    return sum(max(0, min(end, a_end) - max(start, a_start))
               for a_start, a_end in intervals)

# -------------------------
# BED helpers
# -------------------------
def read_bed_subset(bed_file, contig):
    windows = []
    min_start = None
    with open(bed_file) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            chrom, start = fields[0], int(fields[1])
            if chrom != contig:
                continue
            windows.append(fields)
            if min_start is None or start < min_start:
                min_start = start
    return windows, min_start

def index_bed_files(bed_dir, suffix):
    bed_index = {}
    for fname in sorted(os.listdir(bed_dir)):
        if not fname.endswith(f".{suffix}.bed"):
            continue
        base = fname.rsplit(f".{suffix}.bed", 1)[0]
        parts = base.split("_")
        if len(parts) != 2:
            continue
        s1, s2 = parts
        bed_index[frozenset((s1, s2))] = os.path.join(bed_dir, fname)
    return bed_index

# -------------------------
# Main wrapper
# -------------------------
def run_cigar_bed_wrapper(csv_file, bed_dir, bed_suffix, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    bed_index = index_bed_files(bed_dir, bed_suffix)

    with open(csv_file) as csvfh:
        reader = csv.reader(csvfh)
        next(reader)  # skip header row

        for row in reader:
            # sample1, sample1_contig, sample2, sample2_contig, cigar_path
            s1, c1, s2, c2, cigar_file = row

            if not os.path.exists(cigar_file):
                raise FileNotFoundError(cigar_file)

            # Infer ref/query from filename
            # pairwise_cigar_sampleA_sampleB.txt
            fname = os.path.basename(cigar_file)
            parts = fname.replace(".txt", "").split("_")
            if len(parts) != 4:
                raise ValueError(f"Invalid cigar filename: {fname}")

            _, _, ref_sample, query_sample = parts

            if frozenset((ref_sample, query_sample)) != frozenset((s1, s2)):
                raise ValueError(
                    f"Cigar samples {ref_sample},{query_sample} "
                    f"do not match CSV row {s1},{s2}"
                )

            # Map contigs
            ref_contig = c1 if ref_sample == s1 else c2
            query_contig = c2 if query_sample == s2 else c1

            key = frozenset((s1, s2))
            if key not in bed_index:
                raise FileNotFoundError(f"No BED file for {s1}, {s2}")

            bed_file = bed_index[key]

            ref_windows, ref_start = read_bed_subset(bed_file, ref_contig)
            query_windows, query_start = read_bed_subset(bed_file, query_contig)

            if ref_start is None or query_start is None:
                raise ValueError("Missing BED entries for contigs")

            with open(cigar_file) as f:
                cigar = f.readline().strip()

            parsed = parse_cigar(cigar)
            ref_intervals = build_aligned_intervals_from_ops(parsed, ref_start)
            query_intervals = build_aligned_intervals_from_ops(
                reverse_cigar(parsed), query_start
            )

            out_path = os.path.join(
                output_dir, f"{s1}_{s2}.aligned_bases.bed"
            )

            with open(out_path, "w") as out:
                for fields in ref_windows:
                    start, end = int(fields[1]), int(fields[2])
                    out.write("\t".join(fields) + f"\t{count_aligned_bases(ref_intervals, start, end)}\n")
                for fields in query_windows:
                    start, end = int(fields[1]), int(fields[2])
                    out.write("\t".join(fields) + f"\t{count_aligned_bases(query_intervals, start, end)}\n")

            print(f"Processed {s1} vs {s2}")

# -------------------------
# CLI
# -------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Report aligned bases per BED window from pairwise CIGAR string."
    )
    parser.add_argument("-c", "--csv", required=True, help="CSV with columns sample1,sample1_contig,sample2,sample2_contig,cigar path")
    parser.add_argument("-b", "--bed_dir", required=True, help="Directory containing BED files")
    parser.add_argument("-s", "--bed_suffix", required=True, help="BED suffix")
    parser.add_argument("-o", "--out_dir", required=True, help="Output directory")

    args = parser.parse_args()
    run_cigar_bed_wrapper(args.csv, args.bed_dir, args.bed_suffix, args.out_dir)

if __name__ == "__main__":
    main()
