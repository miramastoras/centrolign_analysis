#!/usr/bin/env python3
"""
Unit tests for pairwise_titv_from_cigars.py
"""

import os
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(__file__))
from pairwise_titv_from_cigars import (
    get_mismatch_bases,
    get_near_indel_flags,
    load_fasta,
    parse_cigar,
)

SCRIPT = os.path.join(os.path.dirname(__file__), "pairwise_titv_from_cigars.py")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def write_fasta(path, name, seq):
    with open(path, "w") as fh:
        fh.write(">{}\n{}\n".format(name, seq))


def write_cigar(path, cigar_str):
    with open(path, "w") as fh:
        fh.write(cigar_str)


def run_script(args):
    result = subprocess.run(
        [sys.executable, SCRIPT] + args,
        capture_output=True, text=True
    )
    return result


# ---------------------------------------------------------------------------
# parse_cigar
# ---------------------------------------------------------------------------

class TestParseCigar(unittest.TestCase):

    def test_match_only(self):
        self.assertEqual(parse_cigar("5="), [("=", 5)])

    def test_mismatch(self):
        self.assertEqual(parse_cigar("2=1X3="), [("=", 2), ("X", 1), ("=", 3)])

    def test_indels(self):
        self.assertEqual(parse_cigar("3=1I2=1D1="),
                         [("=", 3), ("I", 1), ("=", 2), ("D", 1), ("=", 1)])

    def test_multi_digit_lengths(self):
        ops = parse_cigar("100=20X5I")
        self.assertEqual(ops, [("=", 100), ("X", 20), ("I", 5)])

    def test_empty_string(self):
        self.assertEqual(parse_cigar(""), [])


# ---------------------------------------------------------------------------
# load_fasta
# ---------------------------------------------------------------------------

class TestLoadFasta(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def test_single_line_seq(self):
        path = os.path.join(self.tmp, "s.fa")
        write_fasta(path, "sample", "ACGT")
        self.assertEqual(load_fasta(path), "ACGT")

    def test_multiline_seq(self):
        path = os.path.join(self.tmp, "s.fa")
        with open(path, "w") as fh:
            fh.write(">sample\nACGT\nACGT\n")
        self.assertEqual(load_fasta(path), "ACGTACGT")

    def test_lowercase_uppercased(self):
        path = os.path.join(self.tmp, "s.fa")
        with open(path, "w") as fh:
            fh.write(">sample\nacgt\n")
        self.assertEqual(load_fasta(path), "ACGT")

    def test_only_first_record_loaded(self):
        path = os.path.join(self.tmp, "s.fa")
        with open(path, "w") as fh:
            fh.write(">rec1\nAAAA\n>rec2\nCCCC\n")
        self.assertEqual(load_fasta(path), "AAAA")


# ---------------------------------------------------------------------------
# get_mismatch_bases
# ---------------------------------------------------------------------------

class TestGetMismatchBases(unittest.TestCase):

    def test_all_matches_no_mismatches(self):
        ops = parse_cigar("5=")
        result = get_mismatch_bases(ops, "AAAAA", "AAAAA")
        self.assertEqual(result, [])

    def test_single_mismatch(self):
        # seq1: AAGAA, seq2: AACAA  →  pos2: G vs C
        ops = parse_cigar("2=1X2=")
        result = get_mismatch_bases(ops, "AAGAA", "AACAA")
        self.assertEqual(result, [("G", "C")])

    def test_multiple_mismatches(self):
        # seq1: AAGAGA, seq2: AAAATA  →  pos2: G/A (Ti), pos4: G/T (Tv)
        ops = parse_cigar("2=1X1=1X1=")
        result = get_mismatch_bases(ops, "AAGAGA", "AAAATA")
        self.assertEqual(result, [("G", "A"), ("G", "T")])

    def test_insertion_before_mismatch_advances_query(self):
        # s1: "AAGAAA" (ref), s2: "AATAAAA" (query, T inserted at qry[2])
        # CIGAR: 2=1I1X3=
        #   2=: ref[0,1]=AA, qry[0,1]=AA
        #   1I: qry[2]='T' consumed, ref stays at 2
        #   1X: ref[2]='G', qry[3]='A'  → Ti
        #   3=: ref[3,4,5]=AAA, qry[4,5,6]=AAA
        ops = parse_cigar("2=1I1X3=")
        result = get_mismatch_bases(ops, "AAGAAA", "AATAAAA")
        self.assertEqual(result, [("G", "A")])

    def test_deletion_before_mismatch_advances_ref(self):
        # s1: "AAAAGAA" (ref, A deleted from query between pos2 and pos3)
        # s2: "AAAGAA" (query)
        # CIGAR: 2=1D1X3=
        #   2=: ref[0,1]=AA, qry[0,1]=AA
        #   1D: ref advances to 3, qry stays at 2
        #   1X: ref[3]='A', qry[2]='G'  → Ti (A/G)
        #   3=: ref[4,5,6]=GAA, qry[3,4,5]=GAA  (wait, need seq to actually match)
        # Let's construct: s1="AAAGAAA", s2="AAAAAA"
        # 2=1D1X3=
        #   2=: AA/AA ✓
        #   1D: ref→3
        #   1X: ref[3]='G', qry[2]='A' → Ti
        #   3=: ref[4,5,6]=AAA, qry[3,4,5]=AAA ✓
        ops = parse_cigar("2=1D1X3=")
        result = get_mismatch_bases(ops, "AAAGAAA", "AAAAAA")
        self.assertEqual(result, [("G", "A")])

    def test_no_ops_returns_empty(self):
        ops = parse_cigar("")
        result = get_mismatch_bases(ops, "AAAA", "AAAA")
        self.assertEqual(result, [])


# ---------------------------------------------------------------------------
# get_near_indel_flags
# ---------------------------------------------------------------------------

class TestGetNearIndelFlags(unittest.TestCase):

    def test_no_mismatches_returns_empty(self):
        ops = parse_cigar("5=")
        self.assertEqual(get_near_indel_flags(ops, 20), [])

    def test_no_indels_all_false(self):
        # One mismatch at ref_pos=3, no indels
        ops = parse_cigar("3=1X3=")
        flags = get_near_indel_flags(ops, 20)
        self.assertEqual(flags, [False])

    def test_snv_adjacent_to_insertion_flagged(self):
        # 2=1I1X3=  insertion at ref_pos=2, mismatch at ref_pos=2  → within window
        ops = parse_cigar("2=1I1X3=")
        flags = get_near_indel_flags(ops, 20)
        self.assertEqual(flags, [True])

    def test_snv_far_from_insertion_not_flagged(self):
        # 2=1I50=1X3=  insertion at ref_pos=2, mismatch at ref_pos=52
        ops = parse_cigar("2=1I50=1X3=")
        flags = get_near_indel_flags(ops, 20)
        self.assertEqual(flags, [False])

    def test_snv_adjacent_to_deletion_flagged(self):
        # 2=1D1X3=  deletion anchors at ref_pos=2, mismatch at ref_pos=3
        ops = parse_cigar("2=1D1X3=")
        flags = get_near_indel_flags(ops, 20)
        self.assertEqual(flags, [True])

    def test_snv_far_from_deletion_not_flagged(self):
        # 2=1D50=1X=  deletion at ref_pos=2, mismatch at ref_pos=53
        ops = parse_cigar("2=1D50=1X1=")
        flags = get_near_indel_flags(ops, 20)
        self.assertEqual(flags, [False])

    def test_window_boundary_exactly_at_limit(self):
        # window=5; insertion at ref_pos=2; mismatch at ref_pos=7  → distance=5 → flagged
        ops = parse_cigar("2=1I4=1X1=")
        flags = get_near_indel_flags(ops, 5)
        self.assertEqual(flags, [True])

    def test_window_boundary_just_outside(self):
        # window=5; insertion at ref_pos=2; mismatch at ref_pos=8  → distance=6 → not flagged
        # 2=1I6=1X1=: after 1I ref stays at 2, after 6= ref=8, mismatch at ref_pos=8
        ops = parse_cigar("2=1I6=1X1=")
        flags = get_near_indel_flags(ops, 5)
        self.assertEqual(flags, [False])

    def test_multiple_snvs_mixed_near_and_far(self):
        # 1=1I1=1X50=1X1=
        # Insertion at ref_pos=1; mismatch1 at ref_pos=2 (near); mismatch2 at ref_pos=52 (far)
        ops = parse_cigar("1=1I1=1X50=1X1=")
        flags = get_near_indel_flags(ops, 20)
        self.assertEqual(flags, [True, False])

    def test_window_zero_flags_only_exact_position(self):
        # window=0: only flag a mismatch at the exact same ref position as an indel.
        # Insertion anchored at ref_pos=2; mismatch at ref_pos=2 → distance=0 → flagged.
        # Note: the main script bypasses get_near_indel_flags entirely when window=0
        # and treats all SNVs as unfiltered (tested in the integration tests).
        ops = parse_cigar("2=1I1X3=")
        flags = get_near_indel_flags(ops, 0)
        self.assertEqual(flags, [True])


# ---------------------------------------------------------------------------
# Ti/Tv classification (via get_mismatch_bases + TRANSITIONS set)
# ---------------------------------------------------------------------------

class TestTiTvClassification(unittest.TestCase):
    """Test that each known transition/transversion pair is classified correctly."""

    def _classify(self, b1, b2):
        from pairwise_titv_from_cigars import TRANSITIONS
        return "Ti" if frozenset([b1, b2]) in TRANSITIONS else "Tv"

    def test_A_G_is_transition(self):
        self.assertEqual(self._classify("A", "G"), "Ti")

    def test_G_A_is_transition(self):
        self.assertEqual(self._classify("G", "A"), "Ti")

    def test_C_T_is_transition(self):
        self.assertEqual(self._classify("C", "T"), "Ti")

    def test_T_C_is_transition(self):
        self.assertEqual(self._classify("T", "C"), "Ti")

    def test_A_C_is_transversion(self):
        self.assertEqual(self._classify("A", "C"), "Tv")

    def test_A_T_is_transversion(self):
        self.assertEqual(self._classify("A", "T"), "Tv")

    def test_G_C_is_transversion(self):
        self.assertEqual(self._classify("G", "C"), "Tv")

    def test_G_T_is_transversion(self):
        self.assertEqual(self._classify("G", "T"), "Tv")

    def test_C_A_is_transversion(self):
        self.assertEqual(self._classify("C", "A"), "Tv")

    def test_T_G_is_transversion(self):
        self.assertEqual(self._classify("T", "G"), "Tv")


# ---------------------------------------------------------------------------
# Integration tests (subprocess)
# ---------------------------------------------------------------------------

class TestIntegration(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.cigar_dir = os.path.join(self.tmp, "cigars")
        self.cigar_prefix = os.path.join(self.cigar_dir, "pairwise_cigar_")
        self.fasta_dir = os.path.join(self.tmp, "fastas")
        os.makedirs(self.cigar_dir)
        os.makedirs(self.fasta_dir)

    def _write_pair(self, s1, s2, seq1, seq2, cigar, chrom="chr1"):
        write_fasta(os.path.join(self.fasta_dir, "{}_{}_hor_array.fasta".format(s1, chrom)),
                    s1, seq1)
        write_fasta(os.path.join(self.fasta_dir, "{}_{}_hor_array.fasta".format(s2, chrom)),
                    s2, seq2)
        write_cigar("{}{}_{}.txt".format(self.cigar_prefix, s1, s2), cigar)

    def _parse_tsv(self, stdout):
        rows = {}
        for line in stdout.strip().splitlines():
            if line.startswith("sample1"):
                continue
            parts = line.split("\t")
            key = (parts[0], parts[1])
            rows[key] = {
                "ti": int(parts[2]),
                "tv": int(parts[3]),
                "ratio": parts[4],
                "filtered": int(parts[5]),
            }
        return rows

    def test_single_transition_A_to_G(self):
        # 2=1X3=, mismatch pos2: ref=G, qry=A  → Ti
        self._write_pair("S1", "S2", "AAGAAA", "AAAAAA", "2=1X3=")
        r = run_script([self.cigar_prefix, self.fasta_dir, "chr1", "--indel_window", "0"])
        rows = self._parse_tsv(r.stdout)
        self.assertEqual(rows[("S1", "S2")]["ti"], 1)
        self.assertEqual(rows[("S1", "S2")]["tv"], 0)
        self.assertEqual(rows[("S1", "S2")]["ratio"], "NA")

    def test_single_transversion_G_to_T(self):
        # 3=1X2=, mismatch pos3: ref=G, qry=T  → Tv
        self._write_pair("S1", "S2", "AAAGAA", "AAATAA", "3=1X2=")
        r = run_script([self.cigar_prefix, self.fasta_dir, "chr1", "--indel_window", "0"])
        rows = self._parse_tsv(r.stdout)
        self.assertEqual(rows[("S1", "S2")]["ti"], 0)
        self.assertEqual(rows[("S1", "S2")]["tv"], 1)

    def test_one_ti_one_tv_ratio_is_1(self):
        # 2=1X1=1X1=
        # pos2: G/A → Ti;  pos4: G/T → Tv
        self._write_pair("S1", "S2", "AAGAGA", "AAAATA", "2=1X1=1X1=")
        r = run_script([self.cigar_prefix, self.fasta_dir, "chr1", "--indel_window", "0"])
        rows = self._parse_tsv(r.stdout)
        self.assertEqual(rows[("S1", "S2")]["ti"], 1)
        self.assertEqual(rows[("S1", "S2")]["tv"], 1)
        self.assertEqual(rows[("S1", "S2")]["ratio"], "1.000000")

    def test_snv_near_insertion_filtered(self):
        # 2=1I1X3=, indel_window=20
        # insertion at ref_pos=2, mismatch at ref_pos=2 → filtered
        self._write_pair("S1", "S2", "AAGAAA", "AATAAAA", "2=1I1X3=")
        r = run_script([self.cigar_prefix, self.fasta_dir, "chr1", "--indel_window", "20"])
        rows = self._parse_tsv(r.stdout)
        self.assertEqual(rows[("S1", "S2")]["ti"], 0)
        self.assertEqual(rows[("S1", "S2")]["tv"], 0)
        self.assertEqual(rows[("S1", "S2")]["filtered"], 1)

    def test_snv_near_insertion_not_filtered_when_window_zero(self):
        # same as above but indel_window=0 disables filtering
        self._write_pair("S1", "S2", "AAGAAA", "AATAAAA", "2=1I1X3=")
        r = run_script([self.cigar_prefix, self.fasta_dir, "chr1", "--indel_window", "0"])
        rows = self._parse_tsv(r.stdout)
        self.assertEqual(rows[("S1", "S2")]["ti"], 1)   # G/A = Ti
        self.assertEqual(rows[("S1", "S2")]["filtered"], 0)

    def test_snv_far_from_insertion_not_filtered(self):
        # insertion at ref_pos=2, then 50 matches, then mismatch at ref_pos=52
        # s1: AA + A*50 + G + A  (53 bases)
        # s2: AA + T (inserted) + A*50 + A + A  (54 bases)
        s1 = "AA" + "A" * 50 + "GA"
        s2 = "AAT" + "A" * 50 + "AA"
        self._write_pair("S1", "S2", s1, s2, "2=1I50=1X1=")
        r = run_script([self.cigar_prefix, self.fasta_dir, "chr1", "--indel_window", "20"])
        rows = self._parse_tsv(r.stdout)
        # mismatch is G/A = Ti, not near indel
        self.assertEqual(rows[("S1", "S2")]["ti"], 1)
        self.assertEqual(rows[("S1", "S2")]["filtered"], 0)

    def test_no_mismatches_outputs_zero_counts(self):
        self._write_pair("S1", "S2", "AAAAA", "AAAAA", "5=")
        r = run_script([self.cigar_prefix, self.fasta_dir, "chr1"])
        rows = self._parse_tsv(r.stdout)
        self.assertEqual(rows[("S1", "S2")]["ti"], 0)
        self.assertEqual(rows[("S1", "S2")]["tv"], 0)
        self.assertIn(rows[("S1", "S2")]["ratio"], ("NA", ""))

    def test_missing_fasta_emits_warning(self):
        # Only write one FASTA; the other is missing
        write_fasta(os.path.join(self.fasta_dir, "S1_chr1_hor_array.fasta"), "S1", "AAGAAA")
        write_cigar("{}S1_S2.txt".format(self.cigar_prefix), "2=1X3=")
        # S2 fasta is missing
        r = run_script([self.cigar_prefix, self.fasta_dir, "chr1"])
        self.assertIn("WARNING", r.stderr)
        self.assertIn("S2", r.stderr)
        self.assertEqual(r.returncode, 0)  # should not crash

    def test_distance_csv_threshold_filters_pairs(self):
        # Two pairs: S1/S2 (distance 0.1, below threshold) and S3/S4 (distance 0.5, above)
        self._write_pair("S1", "S2", "AAGAAA", "AAAAAA", "2=1X3=")
        self._write_pair("S3", "S4", "AACAAA", "AATAAA", "2=1X3=")
        csv_path = os.path.join(self.tmp, "dist.csv")
        with open(csv_path, "w") as fh:
            fh.write("S1,S2,0.1\nS3,S4,0.5\n")
        r = run_script([self.cigar_prefix, self.fasta_dir, "chr1",
                        "--distance_csv", csv_path, "--threshold", "0.2"])
        rows = self._parse_tsv(r.stdout)
        self.assertIn(("S1", "S2"), rows)
        self.assertNotIn(("S3", "S4"), rows)

    def test_distance_csv_includes_both_when_both_below_threshold(self):
        self._write_pair("S1", "S2", "AAGAAA", "AAAAAA", "2=1X3=")
        self._write_pair("S3", "S4", "AACAAA", "AATAAA", "2=1X3=")
        csv_path = os.path.join(self.tmp, "dist.csv")
        with open(csv_path, "w") as fh:
            fh.write("S1,S2,0.1\nS3,S4,0.15\n")
        r = run_script([self.cigar_prefix, self.fasta_dir, "chr1",
                        "--distance_csv", csv_path, "--threshold", "0.2"])
        rows = self._parse_tsv(r.stdout)
        self.assertIn(("S1", "S2"), rows)
        self.assertIn(("S3", "S4"), rows)

    def test_duplicate_reverse_cigar_skipped(self):
        # Both pairwise_cigar_S1_S2.txt and pairwise_cigar_S2_S1.txt present.
        # The second one encountered should be skipped with a warning.
        write_fasta(os.path.join(self.fasta_dir, "S1_chr1_hor_array.fasta"), "S1", "AAGAAA")
        write_fasta(os.path.join(self.fasta_dir, "S2_chr1_hor_array.fasta"), "S2", "AAAAAA")
        write_cigar("{}S1_S2.txt".format(self.cigar_prefix), "2=1X3=")
        write_cigar("{}S2_S1.txt".format(self.cigar_prefix), "2=1X3=")
        r = run_script([self.cigar_prefix, self.fasta_dir, "chr1", "--indel_window", "0"])
        rows = self._parse_tsv(r.stdout)
        # Only one row in output (not two)
        self.assertEqual(len(rows), 1)
        # Warning emitted for the duplicate
        self.assertIn("duplicate", r.stderr)
        self.assertIn("Pairs processed: 1", r.stderr)

    def test_summary_stats_printed_to_stderr(self):
        self._write_pair("S1", "S2", "AAGAAA", "AAAAAA", "2=1X3=")
        r = run_script([self.cigar_prefix, self.fasta_dir, "chr1", "--indel_window", "0"])
        self.assertIn("Pairs processed: 1", r.stderr)
        self.assertIn("Total transitions: 1", r.stderr)
        self.assertIn("Total transversions: 0", r.stderr)

    def test_multiple_mismatches_two_ti_one_tv(self):
        # 1=1X1=1X1=1X1=
        # pos1: G/A Ti,  pos3: C/T Ti,  pos5: G/T Tv
        self._write_pair("S1", "S2", "AGACAGA", "AAATATA", "1=1X1=1X1=1X1=")
        r = run_script([self.cigar_prefix, self.fasta_dir, "chr1", "--indel_window", "0"])
        rows = self._parse_tsv(r.stdout)
        self.assertEqual(rows[("S1", "S2")]["ti"], 2)
        self.assertEqual(rows[("S1", "S2")]["tv"], 1)
        self.assertAlmostEqual(float(rows[("S1", "S2")]["ratio"]), 2.0)


if __name__ == "__main__":
    unittest.main()
