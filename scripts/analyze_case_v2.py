#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 17:06:50 2023

Analyze the accuracy of induced pairwise alignments from an MSA relative to the true
simulated alignment and save the results in a (headerless) table

@author: Jordan

Edited by Mira to pass in "induced_aln_dir" from command line, and output the results file
to a directory specified by command line
"""

import sys
import os
import re
import subprocess
import tempfile

tmp_files = []

def tmp_file_name(suffix = None):
    tmp = os.path.join(tempfile._get_default_tempdir(), next(tempfile._get_candidate_names()))
    if suffix is not None:
        tmp += "." + suffix
    tmp_files.append(tmp)
    return tmp

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

def to_cigar_string(cigar):
    return "".join(str(l) + op for op, l in cigar)

if __name__ == "__main__":

    if len(sys.argv) != 6:
        print("usage: ./analyze_case.py case_dir tree_dist_exec compare_truth_exec induced_aln_dir outdir", file = sys.stderr)
        sys.exit(1)

    case_dir = sys.argv[1]
    tree_dist_exec = sys.argv[2]
    compare_truth_exec = sys.argv[3]
    induced_aln_dir = sys.argv[4]
    outdir=sys.argv[5]

    match = re.search(r'msa_(chr[0-9]{1,2}|chrX|chrY)_sim_cases_\d+/case_(\d+)', case_dir)

    if match:
        chr = match.group(1)  # e.g., 'chrX'
        case = f"case_{match.group(2)}"  # e.g., 'case_24'
        print("Chromosome:", chr)
        print("Case ID:", case)
    else:
        print("No match found.")

    samples = []
    for fp in os.listdir(case_dir):
        m = re.match("sim_(.*)_identity.txt", fp)
        if m is not None:
            samples.append(m.group(1))


    tree_dists = {}
    tree_fp = os.path.join(case_dir, "tree.txt")
    raw_table = subprocess.check_output([tree_dist_exec, tree_fp])
    if type(raw_table) == bytes:
        raw_table = raw_table.decode("utf-8")

    table = [row.split() for row in raw_table.split("\n") if not row.endswith("distance") and not len(row) == 0]
    for s1, s2, d in table:
        tree_dists[(s1, s2)] = float(d)
        tree_dists[(s2, s1)] = float(d)

    # handy to make sure everything lines up even though i think i won't actually print the header...
    header = ["case", "distance", "truth_matches", "truth_match_rate", "matches", "match_rate",
              "mismatches", "mismatch_rate", "recall", "precision"]
    output_table = []

    fp_samp_regex = "_([^_]+)_([^_]+).txt$"
    for aln_fp in os.listdir(induced_aln_dir):
        m = re.search(fp_samp_regex, aln_fp)
        assert(m is not None)

        s1 = m.group(1)
        s2 = m.group(2)

        truth_aln_fp = os.path.join(case_dir, "sim_" + s1 + "_" + s2 + "_cigar.txt")

        if not os.path.exists(truth_aln_fp):
            # we have to reverse the truth relative to this induced alignment
            rev_truth_aln_fp = os.path.join(case_dir, "sim_" + s2 + "_" + s1 + "_cigar.txt")
            assert(os.path.exists(rev_truth_aln_fp))

            tmp_aln_fp = tmp_file_name()
            with open(rev_truth_aln_fp) as aln_in:
                cigar = parse_cigar(aln_in.read().strip())
                reverse_cigar(cigar)
                with open(tmp_aln_fp, "w") as aln_out:
                    print(to_cigar_string(cigar), file = aln_out)

            truth_aln_fp = tmp_aln_fp

        id1_fp = os.path.join(case_dir, "sim_" + s1 + "_identity.txt")
        id2_fp = os.path.join(case_dir, "sim_" + s2 + "_identity.txt")
        raw_comparison = subprocess.check_output([compare_truth_exec, id1_fp, id2_fp, truth_aln_fp, os.path.join(induced_aln_dir, aln_fp)])
        if type(raw_comparison) == bytes:
            raw_comparison = raw_comparison.decode("utf-8")

        row = [None for i in range(len(header))]

        row[header.index("case")] = case
        row[header.index("distance")] = tree_dists[(s1, s2)]

        # parse the output of the comparison binary
        for line in raw_comparison.strip().split("\n"):

            if line.startswith("truth matches"):
                row[header.index("truth_matches")] = int(line.split(": ")[1])

            elif line.startswith("truth match rate"):
                row[header.index("truth_match_rate")] = float(line.split(": ")[1])

            elif line.startswith("aln matches"):
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

        output_table.append(row)



    with open(os.path.join(outdir, chr+"_" + case +"_aln_summary_table.txt"), "w") as out:
        for row in output_table:
            print("\t".join(str(v) for v in row), file = out)

    for tmp_file in tmp_files:
        os.remove(tmp_file)
