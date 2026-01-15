#!/usr/bin/env python3
import argparse
from itertools import combinations
import os
import re

# -------------------------
# Helper functions
# -------------------------

def load_bed(file_path):
    """Load BED file into a list of SV dicts."""
    svs = []
    with open(file_path) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.strip().split()
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            size = end - start
            svs.append({"chrom": chrom, "start": start, "end": end, "size": size})
    return svs

def overlaps(sv1, sv2):
    return max(sv1["start"], sv2["start"]) < min(sv1["end"], sv2["end"])

def breakpoint_distance(sv1, sv2):
    return min(
        abs(sv1["start"] - sv2["start"]),
        abs(sv1["end"] - sv2["end"]),
        abs(sv1["start"] - sv2["end"]),
        abs(sv1["end"] - sv2["start"])
    )

def size_fraction_ok(size1, size2, min_size_frac):
    return min(size1, size2) / max(size1, size2) >= min_size_frac

def overlaps_window(sv, win_start, win_end):
    return max(sv["start"], win_start) < min(sv["end"], win_end)


# -------------------------
# Core matching logic
# -------------------------

def match_sv_windowed(sv_rec, candidates, max_dist, min_size_frac=0.7, debug=False):
    """
    Match one SV against another BED using:
    - expanded window (+/- max_dist)
    - contiguous subsets
    - minimize size difference
    - tie-break by closest breakpoint distance
    - filter by min_size_frac after best subset selection
    """

    sv1 = sv_rec["sv"]

    # Expanded window
    win_start = sv1["start"] - max_dist
    win_end   = sv1["end"] + max_dist

    # Collect candidates overlapping the window
    window_hits = [
        c for c in candidates
        if sv1["chrom"] == c["sv"]["chrom"]
        and overlaps_window(c["sv"], win_start, win_end)
    ]

    if not window_hits:
        if debug:
            sv_rec["reason"] = "no SVs overlap expanded window"
        return []

    # Sort by genomic order
    window_hits.sort(key=lambda c: c["sv"]["start"])

    best_subset = None
    best_size_diff = float("inf")
    best_bp_dist = float("inf")

    n = len(window_hits)

    # Enumerate contiguous subsets
    for i in range(n):
        total_size = 0
        subset_bp_dist = float("inf")

        for j in range(i, n):
            sv2 = window_hits[j]["sv"]
            total_size += sv2["size"]

            subset_bp_dist = min(
                subset_bp_dist,
                breakpoint_distance(sv1, sv2)
            )

            size_diff = abs(total_size - sv1["size"])

            # Primary: size difference
            # Secondary: breakpoint distance
            if (
                size_diff < best_size_diff or
                (size_diff == best_size_diff and subset_bp_dist < best_bp_dist)
            ):
                best_size_diff = size_diff
                best_bp_dist = subset_bp_dist
                best_subset = (i, j, total_size)

    if best_subset is None:
        return []

    i, j, total_size = best_subset

    # Check min_size_frac
    size_ratio = min(sv1["size"], total_size) / max(sv1["size"], total_size)
    if size_ratio < min_size_frac:
        if debug:
            sv_rec["reason"] = (
                f"best subset size fraction failed min_size_frac"
            )
        return []

    # Mark matches
    matched = window_hits[i:j+1]
    sv_rec["matched"] = True
    if debug:
        sv_rec["reason"] = (
            f"windowed-contiguous match"
        )

    for c in matched:
        c["matched"] = True
        if debug:
            c["reason"] = sv_rec["reason"]

    return matched


# -------------------------
# Evaluate BED1 ↔ BED2
# -------------------------

def evaluate_beds(svs1, svs2, max_dist, min_size_frac=0.7, debug=False):
    results1 = [{"sv": sv, "matched": False, "reason": "none"} for sv in svs1]
    results2 = [{"sv": sv, "matched": False, "reason": "none"} for sv in svs2]

    # BED1 -> BED2
    for r1 in results1:
        match_sv_windowed(r1, results2, max_dist, min_size_frac, debug=debug)

    # BED2 -> BED1 (symmetric)
    for r2 in results2:
        match_sv_windowed(r2, results1, max_dist, min_size_frac, debug=debug)

    # Assign final reason
    if debug:
        for r in results1 + results2:
            if not r["matched"] and r["reason"] == "none":
                r["reason"] = "no match"

    return results1, results2

# -------------------------
# Metrics
# -------------------------

def compute_metrics(results1, results2):
    TP1 = sum(r["matched"] for r in results1)
    TP2 = sum(r["matched"] for r in results2)
    P = TP2 / len(results2) if results2 else 0
    R = TP1 / len(results1) if results1 else 0
    F1 = (2 * P * R / (P + R)) if (P + R) > 0 else 0
    HM = F1
    return P, R, F1, HM

# -------------------------
# Output
# -------------------------

def write_results(file_path, results, label, debug=False):
    header = "#chrom\tstart\tend\tsize\tmatched\tlabel"
    if debug:
        header += "\treason"
    with open(file_path, "w") as f:
        f.write(header + "\n")
        for r in results:
            line = f"{r['sv']['chrom']}\t{r['sv']['start']}\t{r['sv']['end']}\t{r['sv']['size']}\t"
            line += f"{'correct' if r['matched'] else 'incorrect'}\t{label}"
            if debug:
                line += f"\t{r.get('reason','')}"
            f.write(line + "\n")

def write_bed9(file_path, results):
    """
    Write BED9 file with RGB coloring:
    green = correct, red = incorrect
    """
    with open(file_path, "w") as f:
        for r in results:
            sv = r["sv"]
            if r["matched"]:
                rgb = "0,255,0"   # green
            else:
                rgb = "255,0,0"   # red

            f.write(
                f"{sv['chrom']}\t"
                f"{sv['start']}\t"
                f"{sv['end']}\t"
                f"0\t"            # name
                f"0\t"            # score
                f"+\t"            # strand
                f"0\t"            # thickStart
                f"0\t"            # thickEnd
                f"{rgb}\n"
            )

# -------------------------
# Main
# -------------------------

def main():
    parser = argparse.ArgumentParser(description="SV comparison with corrected overlap and multi-fragment rescue logic")
    parser.add_argument("bed1", help="BED file 1")
    parser.add_argument("bed2", help="BED file 2")
    parser.add_argument("--max_dist", type=int, default=1000, help="Max breakpoint distance for non-overlapping SVs")
    parser.add_argument("--min_size_frac", type=float, default=0.7, help="Minimum size similarity fraction")
    parser.add_argument("--out_prefix", default="sv_compare", help="Output prefix")
    parser.add_argument("--debug", action="store_true", help="Enable debug mode")
    parser.add_argument("--bed9",action="store_true", help="Output BED9 files with RGB coloring (green=correct, red=incorrect)"
)

    args = parser.parse_args()

    svs1 = load_bed(args.bed1)
    svs2 = load_bed(args.bed2)

    results1, results2 = evaluate_beds(svs1, svs2, args.max_dist, args.min_size_frac, debug=args.debug)

    P, R, F1, HM = compute_metrics(results1, results2)

    summary_path = f"{args.out_prefix}_summary.txt"

    with open(summary_path, "w") as sf:
        sf.write("P\tR\tF1\n")
        sf.write(f"{P:.6f}\t{R:.6f}\t{F1:.6f}\n")

    write_results(f"{args.out_prefix}_BED1_results.txt", results1, "BED1", debug=args.debug)
    write_results(f"{args.out_prefix}_BED2_results.txt", results2, "BED2", debug=args.debug)

    if args.bed9:
        write_bed9(f"{args.out_prefix}_BED1.bed9", results1)
        write_bed9(f"{args.out_prefix}_BED2.bed9", results2)

if __name__ == "__main__":
    main()
