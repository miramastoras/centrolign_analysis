#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
This script 
'''

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

def windowed_alignment_proportion(cigars, win, ref_len):
    """
    Compute windowed alignment proportion over a fixed reference length.

    Parameters
    ----------
    cigars : list of parsed cigars  
        Each cigar is a list of (op, length).
    win : int  
        Window size in reference bases.
    ref_len : int  
        Total reference length (constant for all cigars).

    Returns
    -------
    numpy.ndarray
        Array of windowed average alignment proportions.
    """

    # Allocate coverage array
    cov = np.zeros(ref_len, dtype=np.uint32)

    # Fill coverage using vectorized slices
    for cigar in cigars:
        ref_pos = 0
        for op, length in cigar:
            if op in ("M", "=", "X"):
                # Add 1 to the covered slice
                cov[ref_pos : ref_pos + length] += 1
                ref_pos += length
            elif op == "D":
                # Consumes reference but adds no coverage
                ref_pos += length
            else:
                # I, S, H do NOT consume reference
                continue

    # Convert counts → proportions
    total = len(cigars)
    prop = cov / total

    # Compute windowed averages
    n_windows = (ref_len + win - 1) // win
    windowed = np.zeros(n_windows)

    for i in range(n_windows):
        start = i * win
        end = min((i + 1) * win, ref_len)
        windowed[i] = prop[start:end].mean()

    return windowed


def main():
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


    WINDOW_SIZE = 50  # change as needed
    outdir = "/private/groups/patenlab/mira/windowed_alignment_beds"
    os.makedirs(outdir, exist_ok=True)


    # get list of cigar strings for chr 11 subgroup C 
    cigar_files = [os.path.join(cigar_path, f) for f in os.listdir(cigar_path) if f.startswith("pairwise_cigar")]

    # get dictionary of [ref,query]:filepath for all cigar strings that exist in the directory
    existing_cigar_sample_map = map_to_samples(cigar_files)

    sample_list = smp_lengths['sample'].tolist()

    # Convert sample lengths to dict for fast lookup
    sample_len_dict = dict(zip(smp_lengths["sample"], smp_lengths["length"]))

    proportions=[]

    for sample in sample_list:
        print("Processing sample:", sample)

        ref_len = sample_len_dict[sample]
        pairs = [(sample, other) for other in sample_list if other != sample]
        sample_cigars = []

        for pair in pairs:

            if pair not in existing_cigar_sample_map: # if required pair ordering doesn't exist, we need to reverse the cigar so sample is ref
                rev_pair = (pair[1], pair[0])

                assert rev_pair in existing_cigar_sample_map, f"{rev_pair} not in map!" # if reverse cigar doesn't exit, throw error
                with open(existing_cigar_sample_map[rev_pair], "r") as f:
                    r_cigar = parse_cigar(f.read())
                direct_cigar = reverse_cigar(r_cigar)
            else:
                with open(existing_cigar_sample_map[pair], "r") as f:
                    direct_cigar = parse_cigar(f.read())

            sample_cigars.append(direct_cigar)

        # -----------------------------
        # Compute windowed proportions
        # -----------------------------
        window_props = windowed_alignment_proportion(
            cigars=sample_cigars,
            win=WINDOW_SIZE,
            ref_len=ref_len
        )

        # -----------------------------
        # Write BED file
        # -----------------------------
        bed_path = os.path.join(outdir, f"{sample}.windowed_alignment.bed")

        with open(bed_path, "w") as out:
            for i, prop in enumerate(window_props):
                start = i * WINDOW_SIZE
                end = min((i + 1) * WINDOW_SIZE, ref_len)
                out.write(f"{sample}\t{start}\t{end}\t{prop:.6f}\n")

        print(f"Wrote {bed_path}")


if __name__ == '__main__':
    main()
