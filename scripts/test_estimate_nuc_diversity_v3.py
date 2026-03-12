#!/usr/bin/env python3

"""
Unit tests for estimate_nuc_diversity_v3.py
Run with: pytest test_estimate_nuc_diversity_v3.py -v
"""

import os
import sys
import subprocess
import pytest
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from estimate_nuc_diversity_v3 import parse_cigar, count_aligned, get_near_indel_flags

SCRIPT = os.path.join(os.path.dirname(__file__), "estimate_nuc_diversity_v3.py")


# ---------------------------------------------------------------------------
# Helpers for integration tests
# ---------------------------------------------------------------------------

def run_script(snp_mat, cigar_dir, distance_csv, threshold=None,
               indel_window=None, expect_success=True):
    """Run the script and return (stdout dict, stderr string)."""
    cigar_prefix = os.path.join(cigar_dir, "pairwise_cigar_")
    cmd = [sys.executable, SCRIPT, snp_mat, cigar_prefix, distance_csv]
    if threshold is not None:
        cmd += ["--threshold", str(threshold)]
    if indel_window is not None:
        cmd += ["--indel_window", str(indel_window)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if expect_success:
        assert result.returncode == 0, f"Script failed:\n{result.stderr}"
    out = {}
    for line in result.stdout.splitlines():
        if ":" in line:
            key, _, val = line.partition(":")
            out[key.strip()] = val.strip()
    return out, result.stderr


def write_snp_mat(path, rows):
    n_pos = len(rows[0][1])
    header = "\t".join(["sample"] + [f"pos{i}" for i in range(n_pos)])
    lines = [header]
    for sample, alleles in rows:
        tokens = [sample] + ["?" if a is None else a for a in alleles]
        lines.append("\t".join(tokens))
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def write_distance_csv(path, rows):
    with open(path, "w") as f:
        for s1, s2, d in rows:
            f.write(f"{s1},{s2},{d}\n")


def write_cigar_file(cigar_dir, s1, s2, cigar_str):
    fname = os.path.join(cigar_dir, f"pairwise_cigar_{s1}_{s2}.txt")
    with open(fname, "w") as f:
        f.write(cigar_str + "\n")


# ---------------------------------------------------------------------------
# parse_cigar (carried over from v2)
# ---------------------------------------------------------------------------

class TestParseCigar:
    def test_simple_match(self):
        assert parse_cigar("10=") == [("=", 10)]

    def test_mismatch(self):
        assert parse_cigar("5X") == [("X", 5)]

    def test_mixed_ops(self):
        assert parse_cigar("50=10I40=") == [("=", 50), ("I", 10), ("=", 40)]

    def test_all_cigar_ops(self):
        result = parse_cigar("3H2S5=4X1I2D3N")
        assert result == [("H", 3), ("S", 2), ("=", 5), ("X", 4),
                          ("I", 1), ("D", 2), ("N", 3)]

    def test_empty_string(self):
        assert parse_cigar("") == []


# ---------------------------------------------------------------------------
# count_aligned (carried over from v2)
# ---------------------------------------------------------------------------

class TestCountAligned:
    def test_match_only(self):
        assert count_aligned([("=", 100)]) == 100

    def test_mismatch_counted(self):
        assert count_aligned([("X", 10)]) == 10

    def test_M_counted(self):
        assert count_aligned([("M", 30)]) == 30

    def test_indels_not_counted(self):
        assert count_aligned([("=", 50), ("I", 10), ("=", 40)]) == 90

    def test_soft_hard_clips_not_counted(self):
        assert count_aligned([("H", 5), ("S", 3), ("=", 20)]) == 20

    def test_mixed(self):
        # 30= + 5X + 10I + 2D + 15M = 30+5+15 = 50
        cigar = parse_cigar("30=5X10I2D15M")
        assert count_aligned(cigar) == 50


# ---------------------------------------------------------------------------
# get_near_indel_flags — unit tests
# ---------------------------------------------------------------------------

class TestGetNearIndelFlags:

    def test_no_indels_all_false(self):
        """No I or D in cigar → all mismatches pass through."""
        # 5= 3X 2= — mismatches at ref_pos 5,6,7; no indels
        ops = parse_cigar("5=3X2=")
        flags = get_near_indel_flags(ops, window=20)
        assert list(flags) == [False, False, False]

    def test_no_mismatches_empty_result(self):
        """No X operations → empty result regardless of indels."""
        ops = parse_cigar("5=3D5=3I5=")
        flags = get_near_indel_flags(ops, window=20)
        assert len(flags) == 0

    def test_empty_cigar(self):
        flags = get_near_indel_flags([], window=20)
        assert len(flags) == 0

    def test_x_immediately_adjacent_to_insertion(self):
        """
        Cigar: 5= 1I 1X 5=
          ref_pos after 5=: 5
          1I anchor at 5; ref_pos stays 5
          1X mismatch at ref_pos 5; ref_pos → 6
        Distance |5 - 5| = 0 ≤ 20 → filtered.
        """
        ops = parse_cigar("5=1I1X5=")
        flags = get_near_indel_flags(ops, window=20)
        assert list(flags) == [True]

    def test_x_far_beyond_window_of_insertion(self):
        """
        Cigar: 5= 1I 30= 1X 5=
          1I anchor at ref_pos 5
          1X at ref_pos 36 (5+30+1=36)

        Wait, let me trace:
          5=: ref_pos 0→5
          1I: anchor 5; ref_pos stays 5
          30=: ref_pos 5→35
          1X: mismatch at ref_pos 35
        Distance |35 - 5| = 30 > 20 → not filtered.
        """
        ops = parse_cigar("5=1I30=1X5=")
        flags = get_near_indel_flags(ops, window=20)
        assert list(flags) == [False]

    def test_x_at_exact_window_boundary_is_filtered(self):
        """
        Cigar: 5= 1I 20= 1X 5=
          1I anchor at ref_pos 5
          20=: ref_pos 5→25
          1X mismatch at ref_pos 25
        Distance |25 - 5| = 20 ≤ 20 → filtered (boundary is inclusive).
        """
        ops = parse_cigar("5=1I20=1X5=")
        flags = get_near_indel_flags(ops, window=20)
        assert list(flags) == [True]

    def test_x_one_beyond_window_is_not_filtered(self):
        """
        Cigar: 5= 1I 21= 1X 5=
          1I anchor at ref_pos 5
          21=: ref_pos 5→26
          1X mismatch at ref_pos 26
        Distance |26 - 5| = 21 > 20 → not filtered.
        """
        ops = parse_cigar("5=1I21=1X5=")
        flags = get_near_indel_flags(ops, window=20)
        assert list(flags) == [False]

    def test_deletion_x_near_start_edge(self):
        """
        Deletion anchors both edges. X near start edge → filtered.
        Cigar: 5= 10D 5= 1X 5=
          5=: ref_pos 5
          10D: anchors 5 and 14; ref_pos 5→15
          5=: ref_pos 15→20
          1X: mismatch at ref_pos 20
        Nearest anchor: 14. |20-14|=6 ≤ 20 → filtered.
        """
        ops = parse_cigar("5=10D5=1X5=")
        flags = get_near_indel_flags(ops, window=20)
        assert list(flags) == [True]

    def test_deletion_x_near_end_edge(self):
        """
        X placed right after a long deletion, near the end-edge anchor.
        Cigar: 5= 30D 1= 1X 5=
          5=: ref_pos 5
          30D: anchors 5 and 34; ref_pos 5→35
          1=: ref_pos 35→36
          1X: mismatch at ref_pos 36
        Nearest anchor: 34. |36-34|=2 ≤ 20 → filtered.
        """
        ops = parse_cigar("5=30D1=1X5=")
        flags = get_near_indel_flags(ops, window=20)
        assert list(flags) == [True]

    def test_long_deletion_x_far_from_both_edges(self):
        """
        Very long deletion; X is far from both edges.
        Cigar: 5= 100D 30= 1X 5=
          5=: ref_pos 5
          100D: anchors 5 and 104; ref_pos 5→105
          30=: ref_pos 105→135
          1X: mismatch at ref_pos 135
        |135-5|=130 > 20 and |135-104|=31 > 20 → not filtered.
        """
        ops = parse_cigar("5=100D30=1X5=")
        flags = get_near_indel_flags(ops, window=20)
        assert list(flags) == [False]

    def test_multiple_indels_x_near_any_one_filtered(self):
        """
        Two indels; X near either one should be filtered.
        Cigar: 1X 10= 1I 1= 1X 10= 1D 1= 1X 1=
          1X: mismatch at ref_pos 0; ref_pos→1
          10=: ref_pos 1→11
          1I: anchor 11; ref_pos stays 11
          1=: ref_pos 11→12
          1X: mismatch at ref_pos 12; ref_pos→13
          10=: ref_pos 13→23
          1D: anchors 23 and 23; ref_pos 23→24
          1=: ref_pos 24→25
          1X: mismatch at ref_pos 25; ref_pos→26

        Mismatches: 0, 12, 25
        Anchors: 11, 23, 23

        X at 0:  |0-11|=11 ≤ 20 → True  (near I)
        X at 12: |12-11|=1 ≤ 20 → True  (near I)
        X at 25: |25-23|=2 ≤ 20 → True  (near D)
        """
        ops = parse_cigar("1X10=1I1=1X10=1D1=1X1=")
        flags = get_near_indel_flags(ops, window=20)
        assert list(flags) == [True, True, True]

    def test_mixed_some_near_some_far(self):
        """
        Three mismatches; first and third far from indel, middle one close.
        Cigar: 1X 20= 1I 1= 1X 50= 1X 1=
          1X: mismatch at ref_pos 0; ref_pos→1
          20=: ref_pos 1→21
          1I: anchor 21; ref_pos stays 21
          1=: ref_pos 21→22
          1X: mismatch at ref_pos 22; ref_pos→23
          50=: ref_pos 23→73
          1X: mismatch at ref_pos 73; ref_pos→74

        Mismatches: 0, 22, 73
        Anchors: [21]

        X at 0:  |0-21|=21 > 20  → False
        X at 22: |22-21|=1 ≤ 20  → True
        X at 73: |73-21|=52 > 20 → False
        """
        ops = parse_cigar("1X20=1I1=1X50=1X1=")
        flags = get_near_indel_flags(ops, window=20)
        assert list(flags) == [False, True, False]

    def test_window_zero_no_filtering(self):
        """
        window=0 means only exact ref_pos matches are filtered.
        X immediately following I (same ref_pos): distance=0 → filtered.
        X one position away from I: distance=1 > 0 → not filtered.
        Cigar: 5= 1I 1X 1= 1X 5=
          5=: ref_pos 5
          1I: anchor 5; ref_pos stays 5
          1X: mismatch at ref_pos 5; ref_pos→6  → |5-5|=0 ≤ 0 → True
          1=: ref_pos 6→7
          1X: mismatch at ref_pos 7; ref_pos→8  → |7-5|=2 > 0 → False
        """
        ops = parse_cigar("5=1I1X1=1X5=")
        flags = get_near_indel_flags(ops, window=0)
        assert list(flags) == [True, False]

    def test_custom_window_larger(self):
        """Window of 50 flags a mismatch 40 bp from indel."""
        ops = parse_cigar("5=1I40=1X5=")
        # I anchor at 5; X at 46; |46-5|=41 > 20 but ≤ 50
        flags_20 = get_near_indel_flags(ops, window=20)
        flags_50 = get_near_indel_flags(ops, window=50)
        assert list(flags_20) == [False]
        assert list(flags_50) == [True]

    def test_multi_base_x_op(self):
        """
        A single 3X operation produces 3 separate mismatch positions.
        Cigar: 5= 3X 1I 5=
          5=: ref_pos 5
          3X: mismatches at ref_pos 5,6,7; ref_pos→8
          1I: anchor 8; ref_pos stays 8
          5=: ref_pos 8→13

        Distances: |5-8|=3, |6-8|=2, |7-8|=1 — all ≤ 20 → all True
        """
        ops = parse_cigar("5=3X1I5=")
        flags = get_near_indel_flags(ops, window=20)
        assert list(flags) == [True, True, True]


# ---------------------------------------------------------------------------
# Integration: indel filtering changes Ti/Tv but not total_diffs
# ---------------------------------------------------------------------------

class TestIndelFilteringIntegration:

    def test_snv_near_insertion_excluded_from_titv(self, tmp_path):
        """
        Two SNVs: pos0 (Ti, first in cigar, NOT near indel) and
        pos1 (Ti, second in cigar, near indel → filtered).
        With window=20: only pos0 Ti survives → no Tv → Ti/Tv line absent.
        total_diffs is still 2 (unfiltered).

        Cigar: 1X 20= 1I 1= 1X 50= — mismatches at ref_pos 0 and 22.
        Indel at 21. |0-21|=21>20 (keep), |22-21|=1≤20 (filter).
        """
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "1X20=1I1=1X50=")

        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A", "A"]), ("S2", ["G", "G"])])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])

        out, stderr = run_script(snp_mat, cigar_dir, dist_csv, threshold=0.2,
                                 indel_window=20)
        assert int(out["total nucleotide differences"]) == 2   # unfiltered
        assert "Ti/Tv ratio (differences)" not in out          # all Ti, no Tv → not printed
        assert "SNVs filtered (within 20 ref bp of indel): 1" in stderr

    def test_snv_far_from_insertion_included_in_titv(self, tmp_path):
        """
        Single SNV 30 bp from the only indel (window=20) → not filtered.
        Cigar: 5= 1I 30= 1X 5= — X at ref_pos 36, I anchor at 5.
        |36-5|=31>20 → not filtered.
        """
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "5=1I30=1X5=")

        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A"]), ("S2", ["G"])])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])

        out, stderr = run_script(snp_mat, cigar_dir, dist_csv, threshold=0.2,
                                 indel_window=20)
        assert int(out["total nucleotide differences"]) == 1
        assert "SNVs filtered (within 20 ref bp of indel): 0" in stderr

    def test_indel_window_0_disables_filtering(self, tmp_path):
        """
        --indel_window 0 disables filtering entirely; even SNVs immediately
        adjacent to an indel are included in Ti/Tv.
        Cigar: 1X 1I 1X 5= — X at ref_pos 0 and 1, I anchor at 1.
        With window=0: X at 1 (|1-1|=0 filtered), X at 0 (|0-1|=1 not filtered).
        BUT --indel_window 0 disables the filter entirely in the script
        (the `if args.indel_window > 0` guard), so both are kept.
        """
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "1X1I1X5=")

        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A", "A"]), ("S2", ["G", "C"])])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])

        out, stderr = run_script(snp_mat, cigar_dir, dist_csv, threshold=0.2,
                                 indel_window=0)
        assert int(out["total nucleotide differences"]) == 2
        # Both SNVs included: 1 Ti (A→G) and 1 Tv (A→C) → ratio = 1.0
        assert float(out["Ti/Tv ratio (differences)"]) == pytest.approx(1.0)
        assert "SNVs filtered (within 0 ref bp of indel): 0" in stderr

    def test_filtering_changes_titv_ratio(self, tmp_path):
        """
        Three SNVs: two Ti (first and third in cigar) and one Ti (second, near indel).
        The second SNV is near an indel and gets filtered.

        Cigar: 1X 20= 1I 1= 1X 50= 1X 1=
          X at ref_pos: 0, 22, 73
          I anchor: 21
          |0-21|=21>20 (keep), |22-21|=1≤20 (filter), |73-21|=52>20 (keep)

        SNP matrix: S1=["A","A","A"], S2=["G","G","C"]
          pos0: A→G Ti (kept)
          pos1: A→G Ti (filtered)
          pos2: A→C Tv (kept)

        With window=20:  1 Ti + 1 Tv → ratio = 1.0; filtered = 1
        With window=0 (disabled): 2 Ti + 1 Tv → ratio = 2.0; filtered = 0
        """
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "1X20=1I1=1X50=1X1=")

        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [
            ("S1", ["A", "A", "A"]),
            ("S2", ["G", "G", "C"]),
        ])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])

        out_filtered, stderr_filtered = run_script(
            snp_mat, cigar_dir, dist_csv, threshold=0.2, indel_window=20)
        assert float(out_filtered["Ti/Tv ratio (differences)"]) == pytest.approx(1.0)
        assert "SNVs filtered (within 20 ref bp of indel): 1" in stderr_filtered

        out_unfiltered, stderr_unfiltered = run_script(
            snp_mat, cigar_dir, dist_csv, threshold=0.2, indel_window=0)
        assert float(out_unfiltered["Ti/Tv ratio (differences)"]) == pytest.approx(2.0)
        assert "SNVs filtered (within 0 ref bp of indel): 0" in stderr_unfiltered

    def test_total_diffs_unaffected_by_indel_filter(self, tmp_path):
        """
        total nucleotide differences and nucleotide diversity are computed
        without the indel filter regardless of --indel_window.
        """
        cigar_dir = str(tmp_path)
        # X immediately adjacent to I → would be filtered from Ti/Tv
        write_cigar_file(cigar_dir, "S1", "S2", "100=1I1X100=")

        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A"]), ("S2", ["G"])])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])

        # With strong filtering: SNV is filtered from Ti/Tv
        out_20, _ = run_script(snp_mat, cigar_dir, dist_csv, threshold=0.2,
                               indel_window=20)
        # Without filtering
        out_0, _ = run_script(snp_mat, cigar_dir, dist_csv, threshold=0.2,
                              indel_window=0)

        # total_diffs unchanged regardless of filter
        assert int(out_20["total nucleotide differences"]) == 1
        assert int(out_0["total nucleotide differences"]) == 1
        assert float(out_20["nucleotide diversity"]) == pytest.approx(
            float(out_0["nucleotide diversity"]))

    def test_no_cigar_snvs_included_without_filter(self, tmp_path):
        """
        When no cigar file is found, SNVs still count toward total_diffs and
        Ti/Tv (no indel filtering is possible without a cigar).
        """
        cigar_dir = str(tmp_path)
        # No cigar file written

        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A", "A"]), ("S2", ["G", "C"])])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])

        out, _ = run_script(snp_mat, cigar_dir, dist_csv, threshold=0.2,
                            indel_window=20)
        # total_diffs counts the 2 differences even without cigar
        assert int(out["total nucleotide differences"]) == 2
        # Ti/Tv still computed (1 Ti A→G, 1 Tv A→C → ratio 1.0)
        assert float(out["Ti/Tv ratio (differences)"]) == pytest.approx(1.0)

    def test_snv_near_deletion_filtered(self, tmp_path):
        """
        Deletion; X near end edge of deletion is filtered.
        Cigar: 5= 30D 1= 1X 5=
          30D anchors: 5 and 34; ref_pos after D: 35
          1=: ref_pos 35→36
          1X: mismatch at ref_pos 36
        |36-34|=2 ≤ 20 → filtered.
        """
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "5=30D1=1X5=")

        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A"]), ("S2", ["G"])])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])

        out, stderr = run_script(snp_mat, cigar_dir, dist_csv, threshold=0.2,
                                 indel_window=20)
        assert int(out["total nucleotide differences"]) == 1   # unfiltered
        assert "Ti/Tv ratio (differences)" not in out          # 0 Ti after filter
        assert "SNVs filtered (within 20 ref bp of indel): 1" in stderr

    def test_default_indel_window_is_20(self, tmp_path):
        """
        Without --indel_window, default is 20. X at distance 20 is filtered;
        X at distance 21 is not.
        """
        # X at distance 20 from I: cigar 5= 1I 20= 1X 5=
        #   I anchor 5; X at ref_pos 26; |26-5|=21...
        # Let me recalculate: 5= -> ref_pos 5; 1I anchor 5; 20= -> ref_pos 25; 1X at 25
        # |25-5|=20 exactly → filtered
        cigar_dir_20 = str(tmp_path / "c20")
        cigar_dir_20_obj = tmp_path / "c20"
        cigar_dir_20_obj.mkdir()

        write_cigar_file(str(cigar_dir_20_obj), "S1", "S2", "5=1I20=1X5=")

        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A"]), ("S2", ["G"])])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])

        out, stderr = run_script(snp_mat, str(cigar_dir_20_obj), dist_csv, threshold=0.2)
        # Default window=20; |25-5|=20 ≤ 20 → filtered
        assert "SNVs filtered (within 20 ref bp of indel): 1" in stderr
        assert "Ti/Tv ratio (differences)" not in out

    def test_per_position_event_flags_respect_indel_filter(self, tmp_path):
        """
        Ti/Tv event counts (per-position flags) also reflect the filter.
        Two pairs both have a Ti at pos0. With filtering, pos0 Ti is filtered
        for one pair (near-indel) but kept for the other (far).
        Result: pos0 is still a Ti event (from the unfiltered pair).

        Pair S1-S2: cigar 1X 1I 100= → X at ref_pos 0, I anchor at 1.
                    |0-1|=1 ≤ 20 → filtered. Ti at pos0 excluded.
        Pair S1-S3: cigar 100= 1I 1= 1X 5= → X at ref_pos 101.
                    I anchor at 100. |101-100|=1 ≤ 20 → filtered. Ti at pos0 excluded.

        Hmm, both filtered → no Ti events. Let me instead make one pair
        have its X near an indel and the other not.

        Pair S1-S2: cigar 1X 1I 100=  → X at ref_pos 0, I anchor 1: filtered
        Pair S1-S3: cigar 1X 100=     → X at ref_pos 0, no indels: NOT filtered

        Both pairs differ at pos0 (A→G Ti).
        With window=20: S1-S2 contribution filtered, S1-S3 not filtered.
        → pos0 still gets ti_pos_flags=True from S1-S3.
        Ti/Tv events = 1 Ti event, 0 Tv events → Ti/Tv (events) line absent.
        """
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "1X1I100=")   # X near I → filtered
        write_cigar_file(cigar_dir, "S1", "S3", "1X100=")     # X, no indel → kept

        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [
            ("S1", ["A"]),
            ("S2", ["G"]),
            ("S3", ["G"]),
        ])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [
            ("S1", "S2", 0.1),
            ("S1", "S3", 0.1),
            ("S2", "S3", 0.1),
        ])

        out, _ = run_script(snp_mat, cigar_dir, dist_csv, threshold=0.2,
                            indel_window=20)
        # total_diffs = 2 (S1-S2 + S1-S3); S2-S3 no diff
        assert int(out["total nucleotide differences"]) == 2
        # Ti events: pos0 flagged by S1-S3 (not filtered) → 1 Ti event
        # No Tv events → Ti/Tv (events) line absent
        assert "Ti/Tv ratio (events)" not in out


# ---------------------------------------------------------------------------
# Integration: threshold filtering (carried over from v2)
# ---------------------------------------------------------------------------

class TestThresholdFiltering:
    def test_pairs_below_threshold_included(self, tmp_path):
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "100=")
        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A"]), ("S2", ["G"])])
        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])
        out, _ = run_script(snp_mat, cigar_dir, dist_csv, threshold=0.2)
        assert int(out["total nucleotide differences"]) == 1

    def test_pairs_at_threshold_excluded(self, tmp_path):
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "100=")
        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A"]), ("S2", ["G"])])
        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.2)])  # exactly at threshold
        out, _ = run_script(snp_mat, cigar_dir, dist_csv, threshold=0.2)
        assert int(out["total nucleotide differences"]) == 0

    def test_pairs_above_threshold_excluded(self, tmp_path):
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "100=")
        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A"]), ("S2", ["G"])])
        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.5)])
        out, _ = run_script(snp_mat, cigar_dir, dist_csv, threshold=0.2)
        assert int(out["total nucleotide differences"]) == 0

    def test_default_threshold_is_0_2(self, tmp_path):
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "100=")
        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A"]), ("S2", ["G"])])
        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.19)])
        out, _ = run_script(snp_mat, cigar_dir, dist_csv)
        assert int(out["total nucleotide differences"]) == 1


# ---------------------------------------------------------------------------
# Integration: Ti/Tv classification (carried over from v2)
# ---------------------------------------------------------------------------

class TestTiTvClassification:
    @pytest.mark.parametrize("a1,a2", [
        ("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")
    ])
    def test_transition_pairs(self, tmp_path, a1, a2):
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "100=")
        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", [a1]), ("S2", [a2])])
        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])
        out, _ = run_script(snp_mat, cigar_dir, dist_csv)
        assert int(out["total nucleotide differences"]) == 1
        assert "Ti/Tv ratio (differences)" not in out  # no Tv → not printed

    @pytest.mark.parametrize("a1,a2", [
        ("A", "C"), ("A", "T"), ("G", "C"), ("G", "T"),
        ("C", "A"), ("T", "A"), ("C", "G"), ("T", "G"),
    ])
    def test_transversion_pairs(self, tmp_path, a1, a2):
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "100=")
        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", [a1]), ("S2", [a2])])
        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])
        out, _ = run_script(snp_mat, cigar_dir, dist_csv)
        assert int(out["total nucleotide differences"]) == 1

    def test_ti_tv_ratio(self, tmp_path):
        """2 Ti, 1 Tv → ratio 2.0 (no indels in cigar so no filtering)."""
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "100=")
        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [
            ("S1", ["A", "A", "C"]),
            ("S2", ["G", "C", "T"]),
        ])
        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])
        out, _ = run_script(snp_mat, cigar_dir, dist_csv)
        assert float(out["Ti/Tv ratio (differences)"]) == pytest.approx(2.0)


# ---------------------------------------------------------------------------
# Integration: total aligned and nucleotide diversity (carried over from v2)
# ---------------------------------------------------------------------------

class TestTotalAligned:
    def test_single_pair_cigar(self, tmp_path):
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "50=10I40=")  # aligned = 90
        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A"]), ("S2", ["A"])])
        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])
        out, _ = run_script(snp_mat, cigar_dir, dist_csv)
        assert int(out["total aligned pairs"]) == 90

    def test_multiple_pairs_cigar_summed(self, tmp_path):
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "50=")
        write_cigar_file(cigar_dir, "S2", "S3", "75=")
        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A"]), ("S2", ["A"]), ("S3", ["A"])])
        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1), ("S2", "S3", 0.15)])
        out, _ = run_script(snp_mat, cigar_dir, dist_csv)
        assert int(out["total aligned pairs"]) == 125

    def test_nucleotide_diversity(self, tmp_path):
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "100=")
        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A", "C"]), ("S2", ["G", "C"])])
        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])
        out, _ = run_script(snp_mat, cigar_dir, dist_csv)
        assert float(out["nucleotide diversity"]) == pytest.approx(1.0 / 100)


# ---------------------------------------------------------------------------
# Integration: pair order independence (carried over from v2)
# ---------------------------------------------------------------------------

class TestPairOrder:
    def test_csv_reversed_order_finds_cigar(self, tmp_path):
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "80=")
        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A"]), ("S2", ["G"])])
        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S2", "S1", 0.1)])  # reversed in CSV
        out, _ = run_script(snp_mat, cigar_dir, dist_csv)
        assert int(out["total aligned pairs"]) == 80
        assert int(out["total nucleotide differences"]) == 1


# ---------------------------------------------------------------------------
# Integration: missing data handling (carried over from v2)
# ---------------------------------------------------------------------------

class TestMissingData:
    def test_missing_allele_skipped(self, tmp_path):
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "100=")
        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", [None, "A"]), ("S2", ["G", "G"])])
        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])
        out, _ = run_script(snp_mat, cigar_dir, dist_csv)
        assert int(out["total nucleotide differences"]) == 1

    def test_both_missing_skipped(self, tmp_path):
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "100=")
        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", [None]), ("S2", [None])])
        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])
        out, _ = run_script(snp_mat, cigar_dir, dist_csv)
        assert int(out["total nucleotide differences"]) == 0

    def test_missing_cigar_graceful(self, tmp_path):
        cigar_dir = str(tmp_path)
        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A"]), ("S2", ["G"])])
        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])
        out, _ = run_script(snp_mat, cigar_dir, dist_csv)
        assert int(out["total aligned pairs"]) == 0


# ---------------------------------------------------------------------------
# Integration: event counting (carried over from v2)
# ---------------------------------------------------------------------------

class TestEventCounting:
    def test_same_position_two_pairs_counts_once(self, tmp_path):
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "100=")
        write_cigar_file(cigar_dir, "S1", "S3", "100=")
        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A"]), ("S2", ["G"]), ("S3", ["G"])])
        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1), ("S1", "S3", 0.1)])
        out, _ = run_script(snp_mat, cigar_dir, dist_csv)
        assert int(out["total nucleotide differences"]) == 2
        # 1 Ti event, 0 Tv events → Ti/Tv (events) not printed
        assert "Ti/Tv ratio (events)" not in out

    def test_ti_events_and_tv_events_at_different_positions(self, tmp_path):
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "100=")
        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A", "A"]), ("S2", ["G", "C"])])
        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])
        out, _ = run_script(snp_mat, cigar_dir, dist_csv)
        assert float(out["Ti/Tv ratio (events)"]) == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# Integration: comprehensive scenario (carried over from v2, no indels)
# ---------------------------------------------------------------------------

class TestComprehensive:
    def test_three_samples_no_indels(self, tmp_path):
        """
        3 samples, 2 allowed pairs, cigars with no indels.
        With window=20 and no indels, results should match v2.
        """
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "50=")
        write_cigar_file(cigar_dir, "S2", "S3", "75=")

        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [
            ("S1", ["A",  "C",  None]),
            ("S2", ["G",  "T",  "A" ]),
            ("S3", ["A",  "C",  "T" ]),
        ])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [
            ("S1", "S2", 0.10),
            ("S2", "S3", 0.15),
            ("S1", "S3", 0.50),
        ])

        out, _ = run_script(snp_mat, cigar_dir, dist_csv, threshold=0.2)

        assert int(out["total aligned pairs"]) == 125
        assert int(out["total nucleotide differences"]) == 5
        assert float(out["nucleotide diversity"]) == pytest.approx(5 / 125)
        assert float(out["Ti/Tv ratio (differences)"]) == pytest.approx(4.0)
        assert float(out["Ti/Tv ratio (events)"]) == pytest.approx(2.0)
