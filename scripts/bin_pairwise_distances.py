#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Bin pairwise distances between subset and non-subset samples.")
    parser.add_argument("-s", "--samples", required=True, help="Path to full sample list")
    parser.add_argument("-d", "--distances", required=True, help="CSV with pairwise distances: sample1,sample2,distance")
    parser.add_argument("-u", "--subset", required=True, help="Path to subset sample list")
    parser.add_argument("-o", "--output", required=True,
                    help="Path to write collected distances below 0.9-1.0 bin")
    return parser.parse_args()

def bin_distances(distances, bin_size=0.2, max_val=1.0):
    bins = defaultdict(int)
    bin_edges = [round(i * bin_size, 2) for i in range(int(max_val / bin_size) + 1)]

    for val in distances:
        for i in range(len(bin_edges)-1):
            lower = bin_edges[i]
            upper = bin_edges[i+1]
            # last bin includes the upper bound
            if i == len(bin_edges)-2:
                if lower <= val <= upper:
                    bins[f"[{lower:.1f}-{upper:.1f}]"] += 1
                    break
            else:
                if lower <= val < upper:
                    bins[f"[{lower:.1f}-{upper:.1f})"] += 1
                    break
    return bins

def main():
    args = parse_args()

    # Read full and subset sample lists
    with open(args.samples, 'r') as f:
        all_samples = set(line.strip() for line in f)

    with open(args.subset, 'r') as f:
        subset_samples = set(line.strip() for line in f)

    # Sanity check
    if not subset_samples.issubset(all_samples):
        missing = subset_samples - all_samples
        print(f"[WARNING] These subset samples are not in the full sample list: {', '.join(missing)}")

    non_subset_samples = all_samples - subset_samples

    collected_distances = []

    # Read pairwise distances and collect only subset vs non-subset
    with open(args.distances, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)  # skip header if present
        for row in reader:
            s1, s2, dist = row[0], row[1], float(row[2])
            if (s1 in subset_samples and s2 in non_subset_samples) or \
               (s2 in subset_samples and s1 in non_subset_samples):
                collected_distances.append((s1, s2, dist))

    # Bin the distances
    binned_counts = bin_distances([d[2] for d in collected_distances], bin_size=0.1)

    # Write distances less than 0.9 to file
    # Filter and write distances < 0.9 to output file
    with open(args.output, "w") as out_f:
        out_f.write("sample1,sample2,distance\n")
        for s1, s2, dist in collected_distances:
            if dist < 0.9:
                out_f.write(f"{s1},{s2},{dist:.6f}\n")

    print(f"Wrote filtered distances (< 0.9) to: {args.output}")


    # Print results in bin order
    print("Bin\tCount")
    bin_order = [f"[{i/10:.1f}-{(i+1)/10:.1f}{']' if i == 9 else ')'}" for i in range(10)]
    for b in bin_order:
        print(f"{b}\t{binned_counts.get(b, 0)}")

if __name__ == "__main__":
    main()
