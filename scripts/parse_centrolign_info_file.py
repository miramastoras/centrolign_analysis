#!/usr/bin/env python3
import sys

def load_all_samples(path):
    with open(path) as f:
        return set(line.strip() for line in f if line.strip())

def load_clade_file(path):
    """
    Reads file with format:
       filename    sample1,sample2,sample3
    Returns list of tuples: (filename, set_of_samples)
    """
    clades = []
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            filename = parts[0]
            sample_list = parts[1].split(",")
            clades.append((filename, set(sample_list)))
    return clades

def choose_largest_nonoverlapping_clades(clades):
    """
    Greedy selection of largest clades with NO sample overlap.
    Assumes clade sample sets grow hierarchically (larger sets include smaller ones).
    """
    used_samples = set()
    selected = []

    # Sort clades largest → smallest
    clades_sorted = sorted(clades, key=lambda x: len(x[1]), reverse=True)

    for fname, samples in clades_sorted:
        # Check for overlap with already selected clades
        if not (samples & used_samples):
            selected.append((fname, samples))
            used_samples.update(samples)

    return selected, used_samples

def main():
    if len(sys.argv) != 3:
        print("Usage: python largest_non_overlapping_clades.py <clade_file.txt> <all_samples.txt>")
        sys.exit(1)

    clade_file = sys.argv[1]
    all_samples_file = sys.argv[2]

    # Load input
    all_samples = load_all_samples(all_samples_file)
    clades = load_clade_file(clade_file)

    selected, covered_samples = choose_largest_nonoverlapping_clades(clades)

    # Uncovered samples
    missing = all_samples - covered_samples

    # ---- Output ----
    print("Selected largest non-overlapping clades:\n")
    for fname, samples in selected:
        print(f"{fname}\t{','.join(sorted(samples))}")

    if missing:
        print("\nSamples NOT covered by any clade:")
        print(",".join(sorted(missing)))
    else:
        print("\nAll samples covered.")

if __name__ == "__main__":
    main()
