#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 12:55:28 2023

Analyze the accuracy of direct pairwise alignments of simulated sequences and save the
results in a (headerless) table

@author: Jordan
"""

import sys
import os
import re
import subprocess
import tempfile

if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("usage: ./analyze_pair_case.py case_dir compare_truth_exec", file = sys.stderr)
        sys.exit(1)

    case_dir = sys.argv[1]
    compare_truth_exec = sys.argv[2]

    samples = []
    for fp in os.listdir(case_dir):
        m = re.match("sim_(.*)_identity.txt", fp)
        if m is not None:
            samples.append(m.group(1))

    info = open(os.path.join(case_dir, "sim_info.txt")).read()
    if type(info) == bytes:
        info = info.decode("utf-8")
    generations = int(re.search("generations: (\d+)", info).group(1))

    header = ["case", "aligner", "distance", "truth_matches", "truth_match_rate", "matches", "match_rate",
              "mismatches", "mismatch_rate", "recall", "precision"]
    output_table = []

    truth_fp = os.path.join(case_dir, "sim_seq1_seq2_cigar.txt")
    id1_fp = os.path.join(case_dir, "sim_seq1_identity.txt")
    id2_fp = os.path.join(case_dir, "sim_seq2_identity.txt")


    truth_matches = None
    truth_match_rate = None
    # run the truth against itself to get the truth params
    raw_truth_comparison = subprocess.check_output([compare_truth_exec, id1_fp, id2_fp, truth_fp, truth_fp])
    if type(raw_truth_comparison) == bytes:
        raw_truth_comparison = raw_truth_comparison.decode("utf-8")
    for line in raw_truth_comparison.strip().split("\n"):
        if line.startswith("truth matches"):
            truth_matches = int(line.split(": ")[1])
        elif line.startswith("truth match rate"):
            truth_match_rate = float(line.split(": ")[1])

    assert(truth_matches is not None)
    assert(truth_match_rate is not None)

    for aligner in ["centrolign", "unialigner", "rama"]:

        aln_fp = os.path.join(case_dir, "aln_" + aligner + ".txt")

        row = [case_dir, aligner]
        while len(row) < len(header):
            row.append("NA")

        row[header.index("distance")] = 2 * generations
        row[header.index("truth_matches")] = truth_matches
        row[header.index("truth_match_rate")] = truth_match_rate
        output_table.append(row)

        # check that the alignment completed
        if not os.path.exists(aln_fp):
            continue
        if len(open(aln_fp).read().strip()) == 0:
            continue

        raw_comparison = subprocess.check_output([compare_truth_exec, id1_fp, id2_fp, truth_fp, aln_fp])
        if type(raw_comparison) == bytes:
            raw_comparison = raw_comparison.decode("utf-8")

        # parse the output of the comparison binary
        for line in raw_comparison.strip().split("\n"):

            if line.startswith("aln matches"):
                row[header.index("matches")] = int(line.split(": ")[1])

            elif line.startswith("aln match rate"):
                row[header.index("match_rate")] = float(line.split(": ")[1])

            elif line.startswith("aln mismatches"):
                row[header.index("mismatches")] = int(line.split(": ")[1])

            elif line.startswith("aln mismatch rate"):
                row[header.index("mismatch_rate")] = float(line.split(": ")[1])

            elif line.startswith("aln match completeness"):
                row[header.index("recall")] = float(line.split(": ")[1])

            elif line.startswith("aln match accuracy"):
                row[header.index("precision")] = float(line.split(": ")[1])

    with open(os.path.join(case_dir, "aln_summary_table.txt"), "w") as out:
        for row in output_table:
            print("\t".join(str(v) for v in row), file = out)
