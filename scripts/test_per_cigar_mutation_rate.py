#!/usr/bin/env python3

"""
Unit tests for per_cigar_mutation_rate.py
Run with: pytest test_per_cigar_mutation_rate.py -v
"""

import os
import tempfile
import pytest
import pandas as pd

from per_cigar_mutation_rate import (
    parse_cigar,
    percent_diff,
    passes_adj_filter,
    compute_mutation_counts,
    map_to_samples,
    load_asat_bounds,
)


# ---------------------------------------------------------------------------
# parse_cigar
# ---------------------------------------------------------------------------

class TestParseCigar:
    def test_simple(self):
        assert parse_cigar("10=5X3I2D") == [
            ("=", 10), ("X", 5), ("I", 3), ("D", 2)
        ]

    def test_single_op(self):
        assert parse_cigar("100=") == [("=", 100)]

    def test_all_op_types(self):
        assert parse_cigar("5H10S20=3X4I6D") == [
            ("H", 5), ("S", 10), ("=", 20), ("X", 3), ("I", 4), ("D", 6)
        ]

    def test_large_lengths(self):
        ops = parse_cigar("100000=50000X")
        assert ops == [("=", 100000), ("X", 50000)]

    def test_empty_string(self):
        assert parse_cigar("") == []

    def test_m_op(self):
        """M ops can appear in non-extended CIGAR."""
        assert parse_cigar("10M") == [("M", 10)]


# ---------------------------------------------------------------------------
# percent_diff
# ---------------------------------------------------------------------------

class TestPercentDiff:
    def test_equal_values(self):
        assert percent_diff(10, 10) == 0.0

    def test_normal(self):
        # abs(10-20)/max(10,20) = 10/20 = 0.5
        assert percent_diff(10, 20) == 0.5

    def test_symmetric(self):
        assert percent_diff(20, 10) == percent_diff(10, 20)

    def test_first_zero(self):
        assert percent_diff(0, 10) is False

    def test_second_zero(self):
        assert percent_diff(10, 0) is False

    def test_both_zero(self):
        assert percent_diff(0, 0) is False

    def test_large_diff(self):
        # abs(1-100)/max(1,100) = 99/100 = 0.99
        assert percent_diff(1, 100) == pytest.approx(0.99)

    def test_small_diff(self):
        # abs(49-50)/max(49,50) = 1/50 = 0.02
        assert percent_diff(49, 50) == pytest.approx(0.02)


# ---------------------------------------------------------------------------
# passes_adj_filter
# ---------------------------------------------------------------------------

class TestPassesAdjFilter:
    # --- Case 1: flanked by matches ---
    def test_insertion_flanked_by_matches(self):
        assert passes_adj_filter("I", 10, "=", 50, "=", 50) is True

    def test_deletion_flanked_by_matches(self):
        assert passes_adj_filter("D", 10, "=", 50, "=", 50) is True

    def test_flanked_by_X_and_eq(self):
        assert passes_adj_filter("I", 10, "X", 5, "=", 50) is True

    def test_flanked_by_M(self):
        assert passes_adj_filter("D", 10, "M", 50, "M", 50) is True

    # --- Case 2: adjacent I/D with high diff (> 0.1) ---
    def test_insertion_next_to_deletion_high_diff(self):
        # I=5, next D=50 → percent_diff(5,50) = 45/50 = 0.9
        assert passes_adj_filter("I", 5, "=", 100, "D", 50) is True

    def test_insertion_prev_deletion_high_diff(self):
        # I=5, prev D=50 → percent_diff(5,50) = 0.9
        assert passes_adj_filter("I", 5, "D", 50, "=", 100) is True

    def test_deletion_next_to_insertion_high_diff(self):
        # D=50, next I=5 → percent_diff(50,5) = 0.9
        assert passes_adj_filter("D", 50, "=", 100, "I", 5) is True

    def test_deletion_prev_insertion_high_diff(self):
        # D=50, prev I=5 → percent_diff(50,5) = 0.9
        assert passes_adj_filter("D", 50, "I", 5, "=", 100) is True

    # --- Case 2: adjacent I/D with low diff (<= 0.1) ---
    def test_insertion_next_to_deletion_low_diff(self):
        # I=48, next D=50 → percent_diff(48,50) = 2/50 = 0.04
        assert passes_adj_filter("I", 48, "=", 100, "D", 50) is False

    def test_deletion_next_to_insertion_low_diff(self):
        # D=50, next I=48 → percent_diff(50,48) = 0.04
        assert passes_adj_filter("D", 50, "=", 100, "I", 48) is False

    def test_equal_length_adjacent(self):
        # I=10, next D=10 → percent_diff(10,10) = 0.0
        assert passes_adj_filter("I", 10, "=", 100, "D", 10) is False

    # --- Neither case: should not count ---
    def test_insertion_between_insertions(self):
        # I surrounded by I on both sides
        assert passes_adj_filter("I", 10, "I", 5, "I", 5) is False

    def test_deletion_between_deletions(self):
        assert passes_adj_filter("D", 10, "D", 5, "D", 5) is False

    def test_H_op_not_counted(self):
        # H flanked by matches: passes_adj_filter is only called for I/D
        # but if H were passed, it's in "IHS" branch. However compute_mutation_counts
        # only checks I ops, not H/S. Let's confirm the filter itself:
        # H flanked by matches → Case 1 triggers (prev_op in "MX=", next_op in "MX=")
        assert passes_adj_filter("H", 10, "=", 50, "=", 50) is True
        # But in practice H is never checked by compute_mutation_counts

    def test_insertion_at_cigar_start(self):
        """At start, prev_op defaults to M, so if next_op is also M/X/=, Case 1 triggers."""
        assert passes_adj_filter("I", 10, "M", None, "=", 50) is True

    def test_insertion_at_cigar_end(self):
        """At end, next_op defaults to M."""
        assert passes_adj_filter("I", 10, "=", 50, "M", None) is True


# ---------------------------------------------------------------------------
# compute_mutation_counts
# ---------------------------------------------------------------------------

class TestComputeMutationCounts:
    def test_all_matches(self):
        ops = parse_cigar("100=")
        n_snvs, n_indels, n_svs, aligned = compute_mutation_counts(ops)
        assert n_snvs == 0
        assert n_indels == 0
        assert n_svs == 0
        assert aligned == 100

    def test_snvs_only(self):
        # 10= 5X 10= → 5 SNVs, 25 aligned bases
        ops = parse_cigar("10=5X10=")
        n_snvs, n_indels, n_svs, aligned = compute_mutation_counts(ops)
        assert n_snvs == 5
        assert n_indels == 0
        assert n_svs == 0
        assert aligned == 25

    def test_multiple_snv_blocks(self):
        ops = parse_cigar("10=3X5=2X10=")
        n_snvs, n_indels, n_svs, aligned = compute_mutation_counts(ops)
        assert n_snvs == 5  # 3 + 2
        assert aligned == 30  # 10 + 3 + 5 + 2 + 10

    def test_short_insertion_flanked_by_matches(self):
        # 10= 5I 10= → 1 short indel (flanked by matches, Case 1)
        ops = parse_cigar("10=5I10=")
        n_snvs, n_indels, n_svs, aligned = compute_mutation_counts(ops)
        assert n_snvs == 0
        assert n_indels == 1
        assert n_svs == 0
        assert aligned == 20

    def test_short_deletion_flanked_by_matches(self):
        # 10= 5D 10= → 1 short indel
        ops = parse_cigar("10=5D10=")
        n_snvs, n_indels, n_svs, aligned = compute_mutation_counts(ops)
        assert n_indels == 1
        assert n_svs == 0
        assert aligned == 20

    def test_sv_insertion_flanked_by_matches(self):
        # 10= 100I 10= → 1 SV (> 49bp, flanked by matches)
        ops = parse_cigar("10=100I10=")
        n_snvs, n_indels, n_svs, aligned = compute_mutation_counts(ops)
        assert n_indels == 0
        assert n_svs == 1
        assert aligned == 20

    def test_sv_deletion_flanked_by_matches(self):
        ops = parse_cigar("10=100D10=")
        n_snvs, n_indels, n_svs, aligned = compute_mutation_counts(ops)
        assert n_indels == 0
        assert n_svs == 1
        assert aligned == 20

    def test_boundary_49bp_is_short_indel(self):
        ops = parse_cigar("10=49I10=")
        _, n_indels, n_svs, _ = compute_mutation_counts(ops)
        assert n_indels == 1
        assert n_svs == 0

    def test_boundary_50bp_is_sv(self):
        ops = parse_cigar("10=50I10=")
        _, n_indels, n_svs, _ = compute_mutation_counts(ops)
        assert n_indels == 0
        assert n_svs == 1

    def test_adjacent_id_high_diff_both_counted(self):
        # 10= 5I 50D 10= → I has next_op=D, pd(5,50)=0.9>0.1 → count
        #                    D has prev_op=I, pd(50,5)=0.9>0.1 → count
        ops = parse_cigar("10=5I50D10=")
        _, n_indels, n_svs, aligned = compute_mutation_counts(ops)
        assert n_indels == 1  # 5I is short indel
        assert n_svs == 1     # 50D is SV
        assert aligned == 20

    def test_adjacent_id_low_diff_neither_counted(self):
        # 10= 10I 10D 10= → pd(10,10)=0.0 <= 0.1 → neither counted
        ops = parse_cigar("10=10I10D10=")
        _, n_indels, n_svs, _ = compute_mutation_counts(ops)
        assert n_indels == 0
        assert n_svs == 0

    def test_adjacent_id_at_threshold(self):
        # percent_diff must be STRICTLY > 0.1
        # Need pd(a,b) = 0.1 exactly → abs(a-b)/max(a,b) = 0.1
        # e.g. a=9, b=10 → 1/10 = 0.1 exactly → should NOT count
        ops = parse_cigar("10=9I10D10=")
        _, n_indels, n_svs, _ = compute_mutation_counts(ops)
        assert n_indels == 0
        assert n_svs == 0

    def test_h_and_s_not_counted(self):
        # H and S ops should not be counted as indels/SVs
        ops = parse_cigar("5H10S20=10S5H")
        n_snvs, n_indels, n_svs, aligned = compute_mutation_counts(ops)
        assert n_snvs == 0
        assert n_indels == 0
        assert n_svs == 0
        assert aligned == 20

    def test_mixed_everything(self):
        # 10= 3X 5= 20I 10= 5D 5= 100I 10= 2X 10=
        # SNVs: 3 + 2 = 5
        # aligned: 10+3+5+10+5+10+2+10 = 55
        # 20I flanked by = on both sides → short indel (20 <= 49) → count
        # 5D flanked by = on both sides → short indel → count
        # 100I flanked by = on both sides → SV (100 > 49) → count
        ops = parse_cigar("10=3X5=20I10=5D5=100I10=2X10=")
        n_snvs, n_indels, n_svs, aligned = compute_mutation_counts(ops)
        assert n_snvs == 5
        assert n_indels == 2  # 20I + 5D
        assert n_svs == 1     # 100I
        assert aligned == 55

    def test_insertion_at_start(self):
        # I at position 0: prev_op defaults to "M", next_op is "="
        # Case 1 triggers (M in "MX=", = in "MX=")
        ops = parse_cigar("10I20=")
        _, n_indels, n_svs, aligned = compute_mutation_counts(ops)
        assert n_indels == 1
        assert aligned == 20

    def test_insertion_at_end(self):
        # I at last position: prev_op is "=", next_op defaults to "M"
        ops = parse_cigar("20=10I")
        _, n_indels, n_svs, aligned = compute_mutation_counts(ops)
        assert n_indels == 1
        assert aligned == 20

    def test_empty_cigar(self):
        ops = parse_cigar("")
        n_snvs, n_indels, n_svs, aligned = compute_mutation_counts(ops)
        assert n_snvs == 0
        assert n_indels == 0
        assert n_svs == 0
        assert aligned == 0


# ---------------------------------------------------------------------------
# map_to_samples
# ---------------------------------------------------------------------------

class TestMapToSamples:
    def test_basic(self):
        fps = [
            "/data/cigars/HG002.1_HG003.2.txt",
            "/data/cigars/HG001.1_HG002.2.txt",
        ]
        result = map_to_samples(fps)
        assert ("HG002.1", "HG003.2") in result
        assert ("HG001.1", "HG002.2") in result
        assert result[("HG002.1", "HG003.2")] == fps[0]

    def test_single_file(self):
        fps = ["/path/to/HG002.1_HG003.1.txt"]
        result = map_to_samples(fps)
        assert len(result) == 1
        assert ("HG002.1", "HG003.1") in result


# ---------------------------------------------------------------------------
# load_asat_bounds
# ---------------------------------------------------------------------------

class TestLoadAsatBounds:
    def test_loads_bed_files(self, tmp_path):
        # Create fake bed files
        bed1 = tmp_path / "sampleA.1_asat_arrays.bed"
        bed1.write_text("contig1\t0\t1000\tchr1\n")

        bed2 = tmp_path / "sampleB.2_asat_arrays.bed"
        bed2.write_text("contig2\t500\t2000\tchr1\ncontig3\t100\t800\tchr2\n")

        bounds = load_asat_bounds(str(tmp_path))

        assert bounds[("sampleA.1", "chr1")] == (0, 1000)
        assert bounds[("sampleB.2", "chr1")] == (500, 2000)
        assert bounds[("sampleB.2", "chr2")] == (100, 800)

    def test_empty_dir(self, tmp_path):
        # No bed files → should raise (concat on empty list)
        with pytest.raises(ValueError):
            load_asat_bounds(str(tmp_path))


# ---------------------------------------------------------------------------
# Integration: end-to-end with temp files
# ---------------------------------------------------------------------------

class TestIntegration:
    def test_end_to_end(self, tmp_path):
        """Full pipeline with temp CIGAR files, sample list, and ASAT beds."""
        # Create ASAT bed files
        bed_dir = tmp_path / "asat_beds"
        bed_dir.mkdir()
        (bed_dir / "HG002.1_asat_arrays.bed").write_text("ctg1\t0\t2000\tchr1\n")
        (bed_dir / "HG003.1_asat_arrays.bed").write_text("ctg2\t0\t3000\tchr1\n")

        # Create CIGAR files: 50= 3X 10= 5I 20= 100D 10=
        # Expected: SNVs=3, aligned=93 (50+3+10+20+10), short_indels=1 (5I), SVs=1 (100D)
        cigar_dir = tmp_path / "cigars"
        cigar_dir.mkdir()
        cigar_str = "50=3X10=5I20=100D10="
        (cigar_dir / "cigar_HG002.1_HG003.1.txt").write_text(cigar_str)

        # Create sample list
        sample_file = tmp_path / "samples.txt"
        sample_file.write_text("HG002.1\nHG003.1\n")

        # Output
        output_file = tmp_path / "output.tsv"

        # Run main via subprocess to test CLI
        import subprocess
        result = subprocess.run(
            [
                "python", os.path.join(os.path.dirname(__file__), "per_cigar_mutation_rate.py"),
                "-c", str(cigar_dir / "cigar_"),
                "-s", str(sample_file),
                "-o", str(output_file),
                "-a", str(bed_dir),
                "-chr", "chr1",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, f"stderr: {result.stderr}"

        df = pd.read_csv(output_file, sep="\t")
        assert len(df) == 1

        row = df.iloc[0]
        assert row["n_snvs"] == 3
        assert row["n_short_indels"] == 1
        assert row["n_svs"] == 1
        assert row["aligned_bases"] == 93
        assert row["avg_array_len"] == 2500.0  # (2000+3000)/2

        assert row["n_snvs_per_aligned_base"] == pytest.approx(3 / 93)
        assert row["n_snvs_per_aligned_base_per_avg_len"] == pytest.approx(3 / 93 / 2500)
        assert row["n_short_indels_per_aligned_base"] == pytest.approx(1 / 93)
        assert row["n_short_indels_per_aligned_base_per_avg_len"] == pytest.approx(1 / 93 / 2500)
        assert row["n_svs_per_avg_len"] == pytest.approx(1 / 2500)
