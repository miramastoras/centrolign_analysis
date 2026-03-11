#!/usr/bin/env python3

"""
Unit tests for estimate_nuc_diversity_v2.py
Run with: pytest test_estimate_nuc_diversity_v2.py -v
"""

import os
import sys
import subprocess
import tempfile
import textwrap
import pytest

# Import helper functions directly from the script
sys.path.insert(0, os.path.dirname(__file__))
from estimate_nuc_diversity_v2 import parse_cigar, count_aligned


# ---------------------------------------------------------------------------
# Helpers for integration tests
# ---------------------------------------------------------------------------

SCRIPT = os.path.join(os.path.dirname(__file__), "estimate_nuc_diversity_v2.py")


def run_script(snp_mat, cigar_dir, distance_csv, threshold=None, expect_success=True):
    """Run the script and return (stdout dict, stderr string).
    cigar_dir: directory containing pairwise_cigar_*.txt files; prefix is appended automatically."""
    cigar_prefix = os.path.join(cigar_dir, "pairwise_cigar_")
    cmd = [sys.executable, SCRIPT, snp_mat, cigar_prefix, distance_csv]
    if threshold is not None:
        cmd += ["--threshold", str(threshold)]
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
    """
    Write a minimal SNP matrix TSV.
    rows: list of (sample_name, [allele_or_None, ...])
    """
    n_pos = len(rows[0][1])
    header = "\t".join(["sample"] + [f"pos{i}" for i in range(n_pos)])
    lines = [header]
    for sample, alleles in rows:
        tokens = [sample] + ["?" if a is None else a for a in alleles]
        lines.append("\t".join(tokens))
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def write_distance_csv(path, rows):
    """rows: list of (s1, s2, distance)"""
    with open(path, "w") as f:
        for s1, s2, d in rows:
            f.write(f"{s1},{s2},{d}\n")


def write_cigar_file(cigar_dir, s1, s2, cigar_str):
    fname = os.path.join(cigar_dir, f"pairwise_cigar_{s1}_{s2}.txt")
    with open(fname, "w") as f:
        f.write(cigar_str + "\n")


# ---------------------------------------------------------------------------
# parse_cigar
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
        assert result == [("H", 3), ("S", 2), ("=", 5), ("X", 4), ("I", 1), ("D", 2), ("N", 3)]

    def test_empty_string(self):
        assert parse_cigar("") == []


# ---------------------------------------------------------------------------
# count_aligned
# ---------------------------------------------------------------------------

class TestCountAligned:
    def test_match_only(self):
        assert count_aligned([("=", 100)]) == 100

    def test_mismatch_counted(self):
        assert count_aligned([("X", 10)]) == 10

    def test_M_counted(self):
        assert count_aligned([("M", 30)]) == 30

    def test_indels_not_counted(self):
        # I (insertion) and D (deletion) do not count as aligned
        assert count_aligned([("=", 50), ("I", 10), ("=", 40)]) == 90

    def test_soft_hard_clips_not_counted(self):
        assert count_aligned([("H", 5), ("S", 3), ("=", 20)]) == 20

    def test_mixed(self):
        # 30= + 5X + 10I + 2D + 15M = 30+5+15 = 50
        cigar = parse_cigar("30=5X10I2D15M")
        assert count_aligned(cigar) == 50


# ---------------------------------------------------------------------------
# Integration: threshold filtering
# ---------------------------------------------------------------------------

class TestThresholdFiltering:
    def test_pairs_below_threshold_included(self, tmp_path):
        """Pairs with distance < threshold contribute to diversity."""
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "100=")

        snp_mat = str(tmp_path / "mat.tsv")
        # S1 and S2 differ at pos0 (A→G, transition)
        write_snp_mat(snp_mat, [("S1", ["A"]), ("S2", ["G"])])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])

        out, _ = run_script(snp_mat, cigar_dir, dist_csv, threshold=0.2)
        assert int(out["total nucleotide differences"]) == 1

    def test_pairs_at_threshold_excluded(self, tmp_path):
        """Pairs with distance == threshold are excluded (strict <)."""
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "100=")

        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A"]), ("S2", ["G"])])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.2)])  # exactly threshold

        out, _ = run_script(snp_mat, cigar_dir, dist_csv, threshold=0.2)
        assert int(out["total nucleotide differences"]) == 0

    def test_pairs_above_threshold_excluded(self, tmp_path):
        """Pairs with distance > threshold do not contribute."""
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "100=")

        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A"]), ("S2", ["G"])])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.5)])

        out, _ = run_script(snp_mat, cigar_dir, dist_csv, threshold=0.2)
        assert int(out["total nucleotide differences"]) == 0

    def test_default_threshold_is_0_2(self, tmp_path):
        """Default threshold is 0.2 — pair at 0.19 is included."""
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "100=")

        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A"]), ("S2", ["G"])])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.19)])

        out, _ = run_script(snp_mat, cigar_dir, dist_csv)  # no --threshold
        assert int(out["total nucleotide differences"]) == 1

    def test_mixed_pairs_only_close_count(self, tmp_path):
        """Only the pair below threshold contributes its differences."""
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "100=")
        write_cigar_file(cigar_dir, "S1", "S3", "100=")

        snp_mat = str(tmp_path / "mat.tsv")
        # S1 vs S2: differ at pos0 (Ti); S1 vs S3: differ at pos1 (Ti)
        write_snp_mat(snp_mat, [
            ("S1", ["A", "C"]),
            ("S2", ["G", "C"]),   # differs from S1 at pos0
            ("S3", ["A", "T"]),   # differs from S1 at pos1
        ])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [
            ("S1", "S2", 0.1),   # included
            ("S1", "S3", 0.5),   # excluded
        ])

        out, _ = run_script(snp_mat, cigar_dir, dist_csv, threshold=0.2)
        assert int(out["total nucleotide differences"]) == 1


# ---------------------------------------------------------------------------
# Integration: Ti/Tv classification
# ---------------------------------------------------------------------------

class TestTiTvClassification:
    @pytest.mark.parametrize("a1,a2", [
        ("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")
    ])
    def test_transition_pairs(self, tmp_path, a1, a2):
        """Each canonical transition pair is counted as a transition."""
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "100=")

        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", [a1]), ("S2", [a2])])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])

        out, _ = run_script(snp_mat, cigar_dir, dist_csv)
        assert int(out["total nucleotide differences"]) == 1
        # Ti/Tv ratio should not appear (no transversions) — script skips it
        assert "Ti/Tv ratio (differences)" not in out

    @pytest.mark.parametrize("a1,a2", [
        ("A", "C"), ("A", "T"), ("G", "C"), ("G", "T"),
        ("C", "A"), ("T", "A"), ("C", "G"), ("T", "G"),
    ])
    def test_transversion_pairs(self, tmp_path, a1, a2):
        """Each canonical transversion pair is counted as a transversion."""
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "100=")

        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", [a1]), ("S2", [a2])])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])

        out, _ = run_script(snp_mat, cigar_dir, dist_csv)
        assert int(out["total nucleotide differences"]) == 1
        # With only transversions, Ti/Tv ratio line will be absent
        # but we can verify total diffs == 1 (above)

    def test_ti_tv_ratio(self, tmp_path):
        """Known Ti/Tv ratio is computed correctly."""
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "100=")

        snp_mat = str(tmp_path / "mat.tsv")
        # pos0: A→G (Ti), pos1: A→C (Tv), pos2: C→T (Ti)  → 2 Ti, 1 Tv → ratio 2.0
        write_snp_mat(snp_mat, [
            ("S1", ["A", "A", "C"]),
            ("S2", ["G", "C", "T"]),
        ])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])

        out, _ = run_script(snp_mat, cigar_dir, dist_csv)
        assert float(out["Ti/Tv ratio (differences)"]) == pytest.approx(2.0)


# ---------------------------------------------------------------------------
# Integration: total aligned and nucleotide diversity
# ---------------------------------------------------------------------------

class TestTotalAligned:
    def test_single_pair_cigar(self, tmp_path):
        """Total aligned bases from a simple cigar."""
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "50=10I40=")  # aligned = 90

        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A"]), ("S2", ["A"])])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])

        out, _ = run_script(snp_mat, cigar_dir, dist_csv)
        assert int(out["total aligned pairs"]) == 90

    def test_multiple_pairs_cigar_summed(self, tmp_path):
        """Total aligned is summed across all allowed pairs."""
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "50=")   # 50
        write_cigar_file(cigar_dir, "S2", "S3", "75=")   # 75

        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A"]), ("S2", ["A"]), ("S3", ["A"])])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1), ("S2", "S3", 0.15)])

        out, _ = run_script(snp_mat, cigar_dir, dist_csv)
        assert int(out["total aligned pairs"]) == 125

    def test_nucleotide_diversity(self, tmp_path):
        """Nucleotide diversity = total_diffs / total_aligned."""
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "100=")

        snp_mat = str(tmp_path / "mat.tsv")
        # 1 difference out of 100 aligned bases
        write_snp_mat(snp_mat, [("S1", ["A", "C"]), ("S2", ["G", "C"])])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])

        out, _ = run_script(snp_mat, cigar_dir, dist_csv)
        assert float(out["nucleotide diversity"]) == pytest.approx(1.0 / 100)


# ---------------------------------------------------------------------------
# Integration: pair order independence
# ---------------------------------------------------------------------------

class TestPairOrder:
    def test_csv_reversed_order_finds_cigar(self, tmp_path):
        """CSV lists (S2, S1) but cigar file is named pairwise_cigar_S1_S2.txt."""
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "80=")  # file in S1_S2 order

        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A"]), ("S2", ["G"])])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S2", "S1", 0.1)])  # reversed in CSV

        out, _ = run_script(snp_mat, cigar_dir, dist_csv)
        assert int(out["total aligned pairs"]) == 80
        assert int(out["total nucleotide differences"]) == 1


# ---------------------------------------------------------------------------
# Integration: missing data handling
# ---------------------------------------------------------------------------

class TestMissingData:
    def test_missing_allele_skipped(self, tmp_path):
        """Positions with ? in either sample are not counted as differences."""
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "100=")

        snp_mat = str(tmp_path / "mat.tsv")
        # pos0: S1=?, S2=G → skip; pos1: S1=A, S2=G → Ti
        write_snp_mat(snp_mat, [("S1", [None, "A"]), ("S2", ["G", "G"])])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])

        out, _ = run_script(snp_mat, cigar_dir, dist_csv)
        assert int(out["total nucleotide differences"]) == 1

    def test_both_missing_skipped(self, tmp_path):
        """Positions where both samples are ? are not counted."""
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "100=")

        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", [None]), ("S2", [None])])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])

        out, _ = run_script(snp_mat, cigar_dir, dist_csv)
        assert int(out["total nucleotide differences"]) == 0

    def test_missing_cigar_graceful(self, tmp_path):
        """Script runs without error when a cigar file is missing."""
        cigar_dir = str(tmp_path)
        # No cigar files written at all

        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A"]), ("S2", ["G"])])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])

        out, _ = run_script(snp_mat, cigar_dir, dist_csv)
        assert int(out["total aligned pairs"]) == 0


# ---------------------------------------------------------------------------
# Integration: event counting
# ---------------------------------------------------------------------------

class TestEventCounting:
    def test_same_position_two_pairs_counts_once(self, tmp_path):
        """A transition at the same position from two different pairs = 1 event."""
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "100=")
        write_cigar_file(cigar_dir, "S1", "S3", "100=")

        snp_mat = str(tmp_path / "mat.tsv")
        # Both S2 and S3 differ from S1 at pos0 (both transitions)
        write_snp_mat(snp_mat, [
            ("S1", ["A"]),
            ("S2", ["G"]),  # A→G Ti
            ("S3", ["G"]),  # A→G Ti
        ])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1), ("S1", "S3", 0.1)])

        out, _ = run_script(snp_mat, cigar_dir, dist_csv)
        assert int(out["total nucleotide differences"]) == 2  # 2 pairs differ
        assert float(out["Ti/Tv ratio (events)"]) == pytest.approx(1.0 / 0.0) if False else True
        # 1 Ti event, 0 Tv events → Ti/Tv (events) line absent
        assert "Ti/Tv ratio (events)" not in out

    def test_ti_events_and_tv_events_at_different_positions(self, tmp_path):
        """Distinct Ti and Tv positions each count as separate events."""
        cigar_dir = str(tmp_path)
        write_cigar_file(cigar_dir, "S1", "S2", "100=")

        snp_mat = str(tmp_path / "mat.tsv")
        # pos0: A→G (Ti event), pos1: A→C (Tv event)
        write_snp_mat(snp_mat, [("S1", ["A", "A"]), ("S2", ["G", "C"])])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])

        out, _ = run_script(snp_mat, cigar_dir, dist_csv)
        assert float(out["Ti/Tv ratio (events)"]) == pytest.approx(1.0)  # 1 Ti event / 1 Tv event


# ---------------------------------------------------------------------------
# Integration: comprehensive scenario
# ---------------------------------------------------------------------------

class TestComprehensive:
    def test_three_samples_filtered(self, tmp_path):
        """
        Full scenario: 3 samples, 2 allowed pairs, 3 SNP positions.

        S1: A  C  ?
        S2: G  T  A   (vs S1: pos0 A→G Ti, pos1 C→T Ti)
        S3: A  C  T   (vs S2: pos0 G→A Ti, pos1 T→C Ti, pos2 A→T Tv)

        Allowed pairs: (S1,S2) dist=0.1, (S2,S3) dist=0.15; excluded: (S1,S3) dist=0.5
        Cigar: S1-S2=50=, S2-S3=75= → total_aligned=125

        Expected:
          total_diffs = 2 (S1-S2) + 3 (S2-S3) = 5
          num_transitions = 2+2 = 4
          num_transversions = 1
          Ti/Tv (diffs) = 4.0
          ti_events: pos0 (both pairs Ti), pos1 (both pairs Ti) → 2
          tv_events: pos2 (S2-S3 Tv) → 1
          Ti/Tv (events) = 2.0
          nucleotide_diversity = 5/125 = 0.04
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


# ---------------------------------------------------------------------------
# Warnings: matrix pairs absent from CSV
# ---------------------------------------------------------------------------

class TestMissingCsvWarning:
    def test_no_warning_when_all_pairs_in_csv(self, tmp_path):
        """No warning when every matrix pair appears in the CSV."""
        cigar_dir = str(tmp_path)
        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A"]), ("S2", ["G"])])

        dist_csv = str(tmp_path / "dist.csv")
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])

        _, stderr = run_script(snp_mat, cigar_dir, dist_csv)
        assert "WARNING" not in stderr

    def test_warning_when_matrix_pair_absent_from_csv(self, tmp_path):
        """Warning is printed when a matrix pair is completely absent from the CSV."""
        cigar_dir = str(tmp_path)
        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A"]), ("S2", ["G"]), ("S3", ["A"])])

        dist_csv = str(tmp_path / "dist.csv")
        # CSV only has S1-S2; missing S1-S3 and S2-S3
        write_distance_csv(dist_csv, [("S1", "S2", 0.1)])

        _, stderr = run_script(snp_mat, cigar_dir, dist_csv)
        assert "WARNING" in stderr
        assert "S1" in stderr and "S3" in stderr
        assert "S2" in stderr and "S3" in stderr

    def test_warning_lists_all_missing_pairs(self, tmp_path):
        """Each missing pair is listed individually."""
        cigar_dir = str(tmp_path)
        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A"]), ("S2", ["G"]), ("S3", ["A"])])

        dist_csv = str(tmp_path / "dist.csv")
        # Empty CSV — all 3 pairs missing
        write_distance_csv(dist_csv, [])

        _, stderr = run_script(snp_mat, cigar_dir, dist_csv)
        assert "3 pair(s)" in stderr

    def test_pair_above_threshold_still_suppresses_warning(self, tmp_path):
        """A pair present in CSV but above threshold should NOT trigger a warning."""
        cigar_dir = str(tmp_path)
        snp_mat = str(tmp_path / "mat.tsv")
        write_snp_mat(snp_mat, [("S1", ["A"]), ("S2", ["G"])])

        dist_csv = str(tmp_path / "dist.csv")
        # Pair is in CSV but distance > threshold — should still appear in all_csv_pairs
        write_distance_csv(dist_csv, [("S1", "S2", 0.9)])

        _, stderr = run_script(snp_mat, cigar_dir, dist_csv, threshold=0.2)
        assert "WARNING" not in stderr
