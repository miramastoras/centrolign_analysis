#!/usr/bin/env python3
import argparse
from itertools import combinations

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

# -------------------------
# Core matching logic
# -------------------------

def match_sv(sv, candidates, max_dist, min_size_frac, debug=False):
    """
    Match a single SV against a list of candidate SVs.
    Returns a list of candidate SVs that are matched.
    Updates reason if debug=True.
    """

    # ---- FIX: sv is the wrapper, extract the inner SV ----
    sv_rec = sv
    sv = sv_rec["sv"]

    # Step 1: overlapping SVs
    overlapping = [c for c in candidates if sv["chrom"] == c["sv"]["chrom"] and overlaps(sv, c["sv"])]
    print(overlapping)

    if len(overlapping) == 1:
        c = overlapping[0]
        if size_fraction_ok(sv["size"], c["sv"]["size"], min_size_frac):
            c["matched"] = True
            sv_rec["matched"] = True
            if debug:
                sv_rec["reason"] = "single overlap, size fraction OK"
                c["reason"] = "single overlap, size fraction OK"
            return [c]
        else:
            if debug:
                sv_rec["reason"] = f"single overlap, size fraction FAILED: SV size= {sv["size"]}, candidate size = {c["sv"]["size"]} "
            return []

    elif len(overlapping) > 1:
        # Multiple overlapping SVs: find contiguous subsets with best total size
        overlapping.sort(key=lambda x: x["sv"]["start"])
        n = len(overlapping)
        best_dist = float("inf")
        best_subset = None
        best_size = None
        # Consider all contiguous subsets
        for i in range(n):
            total_size = 0
            for j in range(i, n):
                total_size += overlapping[j]["sv"]["size"]
                dist = abs(total_size - sv["size"])
                if dist < best_dist:
                    best_dist = dist
                    best_subset = (i, j)
                    best_size = total_size
        # Check size fraction
        if size_fraction_ok(sv["size"], best_size, min_size_frac):
            i, j = best_subset
            matched_candidates = overlapping[i:j+1]
            sv_rec["matched"] = True
            if debug:
                sv_rec["reason"] = f"multi-overlap N-way rescue, subset {i}-{j}"
            for c in matched_candidates:
                c["matched"] = True
                if debug:
                    c["reason"] = f"multi-overlap N-way rescue, subset {i}-{j}"
            return matched_candidates
        else:
            if debug:
                sv_rec["reason"] = f"multi-overlap, size fraction FAILED: SV size = {sv["size"]}, best size = {best_size}"
            return []

    else:
        # No overlaps: check SVs within max_dist
        nearby = [c for c in candidates if sv["chrom"] == c["sv"]["chrom"] and breakpoint_distance(sv, c["sv"]) <= max_dist]
        matched_candidates = []
        for c in nearby:
            if size_fraction_ok(sv["size"], c["sv"]["size"], min_size_frac):
                sv_rec["matched"] = True
                c["matched"] = True
                matched_candidates.append(c)
                if debug:
                    sv_rec["reason"] = f"no-overlap but within {max_dist}bp, size fraction OK"
                    c["reason"] = f"no-overlap but within {max_dist}bp, size fraction OK"
        if not matched_candidates and debug:
            sv_rec["reason"] = "no match within {max_dist}bp"
        return matched_candidates

# -------------------------
# Evaluate BED1 ↔ BED2
# -------------------------

def evaluate_beds(svs1, svs2, max_dist, min_size_frac, debug=False):
    results1 = [{"sv": sv, "matched": False, "reason": "none"} for sv in svs1]
    results2 = [{"sv": sv, "matched": False, "reason": "none"} for sv in svs2]

    # BED1 -> BED2
    for r1 in results1:
        match_sv(r1, results2, max_dist, min_size_frac, debug=debug)
    # BED2 -> BED1 (symmetric)
    for r2 in results2:
        match_sv(r2, results1, max_dist, min_size_frac, debug=debug)

    # Assign "no match" reason if unmatched
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
    print(f"Precision: {P:.3f}")
    print(f"Recall:    {R:.3f}")
    print(f"F1:        {F1:.3f}")
    print(f"Harmonic:  {HM:.3f}")

    write_results(f"{args.out_prefix}_BED1_results.txt", results1, "BED1", debug=args.debug)
    write_results(f"{args.out_prefix}_BED2_results.txt", results2, "BED2", debug=args.debug)
    if args.bed9:
        write_bed9(f"{args.out_prefix}_BED1.bed9", results1)
        write_bed9(f"{args.out_prefix}_BED2.bed9", results2)

if __name__ == "__main__":
    main()
