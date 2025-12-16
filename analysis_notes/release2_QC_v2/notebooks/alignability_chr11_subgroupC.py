#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## import statements
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import re
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

cigar_path="/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr11/subgroup_C/induced_pairwise_cigars/"
fai="/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr11/combine_final_subgroups/chr11.subgroup_C.fasta.fai"

# Parse fai file into a dataframe of sample, asat length
smp_lengths = pd.read_csv(
    fai,
    sep="\t",
    header=None,
    usecols=[0, 1],       # only contig and length
    names=["sample", "length"]
)

# This function will take in the length of the array, and current base position, and return the new position normalized between 0 and 1 
def normalize(position, length):
    if length == 1:
        return 0.0  # avoid divide by zero
    return (position - 1) / (length - 1)

## Functions to parse cigar, and reverse cigar
def parse_cigar(cigar):
    parsed = []
    for m in re.finditer("(\d+)([HSMIDX=])", cigar):
        parsed.append((m.group(2), int(m.group(1))))
    return parsed

def reverse_cigar(parsed):
    for i in range(len(parsed)):
        if parsed[i][0] == "I":
            parsed[i] = ("D", parsed[i][1])
        elif parsed[i][0] == "D":
            parsed[i] = ("I", parsed[i][1])
    return parsed

def map_to_samples(fps):
    suff_regex = "([a-zA-Z0-9]+.[0-9])_([a-zA-Z0-9]+.[0-9]).txt$"

    sample_map = {}
    for fp in fps:
        m = re.search(suff_regex, fp)
        assert(m is not None)
        key = tuple([m.group(1), m.group(2)])
        sample_map[key] = fp
    
    return sample_map

## Create list of proportion of aligned coverage at every base in reference cigar string

def per_base_alignment_proportion(cigars):
    '''
    Takes in a nested list of cigar strings all aligned to a single reference 
    Returns a list of proportion aligned at every base, where index = base position in array
    '''
    # Build coverage dictionary mapping ref_pos → count of aligned M/= /X
    coverage = defaultdict(int)

    for cigar in cigars:
        ref_pos = 0
        
        for op, length in cigar:
            if op in ("M", "=", "X"):
                for _ in range(length):
                    coverage[ref_pos] += 1
                    ref_pos += 1
            elif op == "D":
                # consumes ref but is not counted as aligned
                ref_pos += length
            else:
                # I, S, H do NOT consume ref
                continue

    # Find reference length
    max_ref_pos = max(coverage.keys()) + 1

    total = len(cigars)
    proportions = [coverage[pos] / total for pos in range(max_ref_pos)]
    return proportions


# get list of cigar strings for chr 11 subgroup C 
cigar_files = [os.path.join(cigar_path, f) for f in os.listdir(cigar_path) if f.startswith("pairwise_cigar")]

# get dictionary of [ref,query]:filepath for all cigar strings in 
existing_cigar_sample_map = map_to_samples(cigar_files)

sample_list = smp_lengths['sample'].tolist()
proportions=[]

# Store all curves for plotting
all_sample_curves = {}

for sample in sample_list:
    print("Processing sample:", sample)
    
    pairs = [(sample, other) for other in sample_list if other != sample]
    sample_cigars = []

    for pair in pairs:

        if pair not in existing_cigar_sample_map:
            rev_pair = (pair[1], pair[0])

            assert rev_pair in existing_cigar_sample_map, f"{rev_pair} not in map!"
            with open(existing_cigar_sample_map[rev_pair], "r") as f:
                r_cigar = parse_cigar(f.read())
            direct_cigar = reverse_cigar(r_cigar)
        else:
            with open(existing_cigar_sample_map[pair], "r") as f:
                direct_cigar = parse_cigar(f.read())

        sample_cigars.append(direct_cigar)

    proportions = per_base_alignment_proportion(sample_cigars)
    all_sample_curves[sample] = proportions

# -------------------------------------------------------------------
# Plotting: one line per sample
# -------------------------------------------------------------------

plt.figure(figsize=(12, 6))

colors = cm.get_cmap('tab20', len(sample_list))  # 'tab20' is a qualitative colormap

for i, sample in enumerate(sample_list):
    proportions = all_sample_curves[sample]
    length = smp_lengths.loc[smp_lengths["sample"] == sample, "length"].values[0]

    # Normalize the x-axis positions
    xs = [normalize(i+1, length) for i in range(len(proportions))]
    ys = proportions

    # Plot line connecting points
    plt.plot(xs, ys, color=colors(i), linewidth=1, alpha=0.7)
    plt.scatter(xs, ys, color=colors(i), s=10, alpha=0.7)

plt.xlabel("Normalized base position (0–1)")
plt.ylabel("Proportion aligned")
plt.title("Per-base alignment proportion across samples")
plt.tight_layout()

# ---- SAVE THE FIGURE ----
output_file = "/private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/plots/alignment_proportion_plot.png"
plt.savefig(output_file, dpi=300)
print("Saved plot to:", output_file)

plt.show()