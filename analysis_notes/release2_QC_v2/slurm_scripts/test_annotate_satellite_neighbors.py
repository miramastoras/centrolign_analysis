#!/usr/bin/env python3
"""Unit tests for annotate_satellite_neighbors.py"""

import os, sys, tempfile, unittest

# import the module under test
sys.path.insert(0, os.path.dirname(__file__))
from annotate_satellite_neighbors import (
    strip_contig, is_active_hor, is_ct, parse_bed, has_large_flanking_blocks,
    merge_nonct_blocks, ACROCENTRIC
)


# ── helpers ───────────────────────────────────────────────────────────────────

def make_bed(lines, tmp_dir):
    p = os.path.join(tmp_dir, 'test.bed')
    with open(p, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    return p

def make_tsv(lines, tmp_dir, name='HG00097_hap1_chr13_intersect.tsv'):
    p = os.path.join(tmp_dir, name)
    with open(p, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    return p

def gene_row(contig, start, end):
    """Minimal TSV row matching the real 11-col format."""
    return f"{contig}\tCAT\tgene\t{start}\t{end}\t.\t+\t.\tattributes\t{contig}\t100"


# ── strip_contig ──────────────────────────────────────────────────────────────

class TestStripContig(unittest.TestCase):

    def test_plain_name_unchanged(self):
        self.assertEqual(strip_contig('chr13'), 'chr13')

    def test_triple_hash_strips_to_last(self):
        self.assertEqual(strip_contig('HG00097#hap1#CM094067.1'), 'CM094067.1')

    def test_double_hash_strips_to_last(self):
        self.assertEqual(strip_contig('sample#CM094067.1'), 'CM094067.1')

    def test_no_hash_returns_as_is(self):
        self.assertEqual(strip_contig('CM094067.1'), 'CM094067.1')


# ── is_active_hor ─────────────────────────────────────────────────────────────

class TestIsActiveHor(unittest.TestCase):

    def test_exact_match(self):
        self.assertTrue(is_active_hor('active_hor'))

    def test_case_insensitive(self):
        self.assertTrue(is_active_hor('Active_HOR'))
        self.assertTrue(is_active_hor('ACTIVE_HOR_1_1'))

    def test_with_suffix(self):
        self.assertTrue(is_active_hor('active_hor_1_1(S1C1H1)'))

    def test_non_active_hor(self):
        self.assertFalse(is_active_hor('hor_1_1'))
        self.assertFalse(is_active_hor('dhor'))
        self.assertFalse(is_active_hor('ct_1_1'))
        self.assertFalse(is_active_hor('HSAT2'))


# ── parse_bed ─────────────────────────────────────────────────────────────────

class TestParseBed(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def test_basic_parsing(self):
        bed = make_bed([
            'chr1\t1000\t200000\tHSAT2\t100\t.\t1000\t200000\t0,0,255',
            'chr1\t300000\t500000\tactive_hor_1\t100\t.\t300000\t500000\t255,0,0',
        ], self.tmp)
        blocks, hor_starts = parse_bed(bed)
        self.assertIn('chr1', blocks)
        self.assertEqual(len(blocks['chr1']), 2)
        self.assertIn('chr1', hor_starts)
        self.assertEqual(hor_starts['chr1'], 300000)

    def test_track_line_skipped(self):
        bed = make_bed([
            'track name=censat',
            'chr1\t1000\t200000\tHSAT2\t100\t.\t1000\t200000\t0,0,255',
        ], self.tmp)
        blocks, _ = parse_bed(bed)
        self.assertEqual(len(blocks['chr1']), 1)

    def test_hash_prefix_stripped(self):
        bed = make_bed([
            'HG00097#hap1#CM094067.1\t1000\t200000\tHSAT2\t100\t.\t1000\t200000\t0,0,255',
        ], self.tmp)
        blocks, _ = parse_bed(bed)
        self.assertIn('CM094067.1', blocks)
        self.assertNotIn('HG00097#hap1#CM094067.1', blocks)

    def test_blocks_sorted_by_start(self):
        bed = make_bed([
            'chr1\t500000\t700000\tHSAT3\t100\t.\t500000\t700000\t0,0,255',
            'chr1\t1000\t200000\tHSAT2\t100\t.\t1000\t200000\t0,0,255',
        ], self.tmp)
        blocks, _ = parse_bed(bed)
        starts = [b[0] for b in blocks['chr1']]
        self.assertEqual(starts, sorted(starts))

    def test_hor_start_takes_minimum(self):
        bed = make_bed([
            'chr13\t500000\t700000\tactive_hor_1\t100\t.\t500000\t700000\t255,0,0',
            'chr13\t100000\t300000\tactive_hor_2\t100\t.\t100000\t300000\t255,0,0',
        ], self.tmp)
        _, hor_starts = parse_bed(bed)
        self.assertEqual(hor_starts['chr13'], 100000)

    def test_no_active_hor_absent_from_hor_starts(self):
        bed = make_bed([
            'chr1\t1000\t200000\tHSAT2\t100\t.\t1000\t200000\t0,0,255',
        ], self.tmp)
        _, hor_starts = parse_bed(bed)
        self.assertNotIn('chr1', hor_starts)

    def test_size_computed_correctly(self):
        bed = make_bed([
            'chr1\t0\t150000\tHSAT2\t100\t.\t0\t150000\t0,0,255',
        ], self.tmp)
        blocks, _ = parse_bed(bed)
        self.assertEqual(blocks['chr1'][0][2], 150000)


# ── is_ct ─────────────────────────────────────────────────────────────────────

class TestIsCt(unittest.TestCase):

    def test_ct_matches(self):
        self.assertTrue(is_ct('ct_1_1'))
        self.assertTrue(is_ct('CT_1_2'))
        self.assertTrue(is_ct('ct'))

    def test_non_ct(self):
        self.assertFalse(is_ct('active_hor'))
        self.assertFalse(is_ct('HSAT2'))
        self.assertFalse(is_ct('hor_1_1'))


# ── merge_nonct_blocks ────────────────────────────────────────────────────────

class TestMergeNonctBlocks(unittest.TestCase):

    def test_empty(self):
        self.assertEqual(merge_nonct_blocks([]), [])

    def test_all_ct_unchanged(self):
        blocks = [(0, 100000, 100000, 'ct_1_1'), (200000, 300000, 100000, 'ct_1_2')]
        result = merge_nonct_blocks(blocks)
        self.assertEqual(len(result), 2)
        self.assertEqual(result[0][3], 'ct_1_1')
        self.assertEqual(result[1][3], 'ct_1_2')

    def test_adjacent_nonct_merged(self):
        blocks = [
            (0,      100000, 100000, 'HSAT2'),
            (100000, 200000, 100000, 'HSAT3'),  # touching
        ]
        result = merge_nonct_blocks(blocks)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], (0, 200000, 200000, 'satellite'))

    def test_nonct_separated_by_ct_not_merged(self):
        blocks = [
            (0,      100000, 100000, 'HSAT2'),
            (100000, 200000, 100000, 'ct_1_1'),
            (200000, 300000, 100000, 'HSAT3'),
        ]
        result = merge_nonct_blocks(blocks)
        self.assertEqual(len(result), 3)
        self.assertEqual(result[0][3], 'satellite')
        self.assertEqual(result[1][3], 'ct_1_1')
        self.assertEqual(result[2][3], 'satellite')

    def test_overlapping_nonct_merged(self):
        blocks = [
            (0,     150000, 150000, 'HSAT2'),
            (100000, 250000, 150000, 'HSAT3'),  # overlaps previous
        ]
        result = merge_nonct_blocks(blocks)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], (0, 250000, 250000, 'satellite'))

    def test_gap_between_nonct_not_merged(self):
        blocks = [
            (0,      100000, 100000, 'HSAT2'),
            (200000, 300000, 100000, 'HSAT3'),  # gap: 100k-200k
        ]
        result = merge_nonct_blocks(blocks)
        self.assertEqual(len(result), 2)

    def test_ct_flanked_by_nonct(self):
        blocks = [
            (0,      100000, 100000, 'HSAT2'),
            (100000, 200000, 100000, 'HSAT3'),
            (200000, 300000, 100000, 'ct_1_1'),
            (300000, 400000, 100000, 'mon'),
            (400000, 500000, 100000, 'dhor'),
        ]
        result = merge_nonct_blocks(blocks)
        self.assertEqual(len(result), 3)
        self.assertEqual(result[0], (0, 200000, 200000, 'satellite'))
        self.assertEqual(result[1][3], 'ct_1_1')
        self.assertEqual(result[2], (300000, 500000, 200000, 'satellite'))

    def test_merged_size_correct(self):
        blocks = [
            (0, 80000, 80000, 'HSAT2'),
            (80000, 160000, 80000, 'HSAT3'),
            (160000, 300000, 140000, 'mon'),
        ]
        result = merge_nonct_blocks(blocks)
        self.assertEqual(result[0][2], 300000)  # size = end - start


# ── has_large_flanking_blocks ─────────────────────────────────────────────────

class TestHasLargeFlankingBlocks(unittest.TestCase):

    def _blocks(self, *intervals, name='HSAT2'):
        return sorted([(s, e, e - s, name) for s, e in intervals])

    # All tests below pin min_size=100_000 to test logic independently of the
    # default threshold (which can be changed without breaking these tests).

    def test_both_sides_large(self):
        blocks = self._blocks((0, 200000), (400000, 700000))
        l, r = has_large_flanking_blocks(blocks, 250000, 350000, min_size=100_000)
        self.assertTrue(l)
        self.assertTrue(r)

    def test_left_too_small(self):
        blocks = self._blocks((200000, 250000), (400000, 700000))  # left = 50kb
        l, r = has_large_flanking_blocks(blocks, 250000, 350000, min_size=100_000)
        self.assertFalse(l)
        self.assertTrue(r)

    def test_right_too_small(self):
        blocks = self._blocks((0, 200000), (400000, 450000))  # right = 50kb
        l, r = has_large_flanking_blocks(blocks, 250000, 350000, min_size=100_000)
        self.assertTrue(l)
        self.assertFalse(r)

    def test_no_blocks_on_either_side(self):
        l, r = has_large_flanking_blocks([], 250000, 350000, min_size=100_000)
        self.assertFalse(l)
        self.assertFalse(r)

    def test_no_left_block(self):
        blocks = self._blocks((400000, 700000),)
        l, r = has_large_flanking_blocks(blocks, 250000, 350000, min_size=100_000)
        self.assertFalse(l)
        self.assertTrue(r)

    def test_no_right_block(self):
        blocks = self._blocks((0, 200000),)
        l, r = has_large_flanking_blocks(blocks, 250000, 350000, min_size=100_000)
        self.assertTrue(l)
        self.assertFalse(r)

    def test_block_touching_gene_start_counts_as_left(self):
        blocks = self._blocks((0, 250000),)  # 250kb, ends at gene_start
        l, _ = has_large_flanking_blocks(blocks, 250000, 350000, min_size=100_000)
        self.assertTrue(l)

    def test_block_touching_gene_end_counts_as_right(self):
        blocks = self._blocks((350000, 600000),)  # 250kb, starts at gene_end
        _, r = has_large_flanking_blocks(blocks, 250000, 350000, min_size=100_000)
        self.assertTrue(r)

    def test_overlapping_block_not_counted(self):
        blocks = self._blocks((200000, 300000),)  # overlaps gene 250k-350k
        l, r = has_large_flanking_blocks(blocks, 250000, 350000, min_size=100_000)
        self.assertFalse(l)
        self.assertFalse(r)

    def test_non_nearest_large_left_block_passes(self):
        # nearest left block is small (50kb), but a farther one is large (200kb)
        blocks = self._blocks((0, 200000), (210000, 250000))  # far=200kb, near=40kb
        l, _ = has_large_flanking_blocks(blocks, 250000, 350000, min_size=100_000)
        self.assertTrue(l)

    def test_non_nearest_large_right_block_passes(self):
        # nearest right block is small (50kb), but a farther one is large (200kb)
        blocks = self._blocks((350000, 400000), (450000, 700000))  # near=50kb, far=250kb
        _, r = has_large_flanking_blocks(blocks, 250000, 350000, min_size=100_000)
        self.assertTrue(r)

    def test_large_ct_block_to_left_does_not_count(self):
        # A large CT block (> min_size) to the left must NOT trigger left_ok —
        # only satellite blocks count. This was a real bug: genes just past the
        # right edge of a large CT block were incorrectly getting left_ok=True.
        blocks = [
            (0,      600000, 600000, 'ct'),        # large CT block to the left
            (600000, 650000, 50000,  'satellite'),  # small satellite between CT and gene
            (800000, 850000, 50000,  'satellite'),  # small satellite to the right
        ]
        l, r = has_large_flanking_blocks(blocks, 660000, 780000, min_size=100_000)
        self.assertFalse(l)   # CT block must not count
        self.assertFalse(r)

    def test_all_small_blocks_both_sides_false(self):
        blocks = self._blocks((0, 50000), (200000, 230000),
                               (400000, 450000), (500000, 530000))
        l, r = has_large_flanking_blocks(blocks, 250000, 350000, min_size=100_000)
        self.assertFalse(l)
        self.assertFalse(r)

    def test_custom_min_size(self):
        blocks = self._blocks((0, 80000), (400000, 480000))   # each 80kb
        l, r = has_large_flanking_blocks(blocks, 250000, 350000, min_size=50000)
        self.assertTrue(l)
        self.assertTrue(r)
        l2, r2 = has_large_flanking_blocks(blocks, 250000, 350000, min_size=100000)
        self.assertFalse(l2)
        self.assertFalse(r2)


# ── end-to-end via main() ─────────────────────────────────────────────────────

class TestEndToEnd(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def _run(self, tsv_lines, bed_lines, tsv_name='HG00097_hap1_chr13_intersect.tsv'):
        import annotate_satellite_neighbors as ann
        bed = make_bed(bed_lines, self.tmp)
        tsv = make_tsv(tsv_lines, self.tmp, name=tsv_name)
        out = os.path.join(self.tmp, 'out.tsv')

        # patch argv and call main
        orig = sys.argv[:]
        sys.argv = ['x', tsv, bed, out]
        try:
            ann.main()
        finally:
            sys.argv = orig

        with open(out) as f:
            lines = [line.rstrip('\n').split('\t') for line in f]
        return lines[1:]  # skip header

    def test_between_true_when_majority_ct_and_both_neighbors_large(self):
        bed_lines = [
            'CM094067.1\t0\t200000\tHSAT2\t100\t.\t0\t200000\t0,0,255',           # 200kb > 100kb
            'CM094067.1\t200000\t400000\tct_1_1\t100\t.\t200000\t400000\t224,224,224',
            'CM094067.1\t400000\t700000\tHSAT3\t100\t.\t400000\t700000\t0,0,255',  # 300kb > 100kb
        ]
        tsv_lines = [gene_row('CM094067.1', 250000, 350000)]
        rows = self._run(tsv_lines, bed_lines)
        self.assertEqual(rows[0][-2], 'TRUE')

    def test_between_false_when_majority_ct_and_left_too_small(self):
        bed_lines = [
            'CM094067.1\t200000\t250000\tHSAT2\t100\t.\t200000\t250000\t0,0,255',  # 50kb < 100kb
            'CM094067.1\t250000\t370000\tct_1_1\t100\t.\t250000\t370000\t224,224,224',
            'CM094067.1\t400000\t700000\tHSAT3\t100\t.\t400000\t700000\t0,0,255',  # 300kb > 100kb
        ]
        tsv_lines = [gene_row('CM094067.1', 260000, 350000)]
        rows = self._run(tsv_lines, bed_lines)
        self.assertEqual(rows[0][-2], 'FALSE')

    def test_between_false_when_non_ct_gene_has_no_large_flanking(self):
        # gene overlaps an HSAT2 block but has no >100kb satellite on either side
        bed_lines = [
            'CM094067.1\t200000\t400000\tHSAT2\t100\t.\t200000\t400000\t0,0,255',
        ]
        tsv_lines = [gene_row('CM094067.1', 250000, 350000)]
        rows = self._run(tsv_lines, bed_lines)
        self.assertEqual(rows[0][-2], 'FALSE')

    def test_between_false_when_no_overlap_at_all(self):
        # gene has no overlapping censat block and no large flanking → FALSE
        bed_lines = [
            'CM094067.1\t0\t100000\tHSAT2\t100\t.\t0\t100000\t0,0,255',
        ]
        tsv_lines = [gene_row('CM094067.1', 500000, 600000)]
        rows = self._run(tsv_lines, bed_lines)
        self.assertEqual(rows[0][-2], 'FALSE')

    def test_between_false_when_right_missing_and_majority_ct(self):
        bed_lines = [
            'CM094067.1\t0\t200000\tHSAT2\t100\t.\t0\t200000\t0,0,255',
            'CM094067.1\t200000\t400000\tct_1_1\t100\t.\t200000\t400000\t224,224,224',
        ]
        tsv_lines = [gene_row('CM094067.1', 250000, 350000)]
        rows = self._run(tsv_lines, bed_lines)
        self.assertEqual(rows[0][-2], 'FALSE')

    def test_acrocentric_arm_true(self):
        # chr13 gene entirely left of active_hor start
        bed_lines = [
            'CM094067.1\t0\t200000\tHSAT2\t100\t.\t0\t200000\t0,0,255',
            'CM094067.1\t500000\t700000\tactive_hor_1\t100\t.\t500000\t700000\t255,0,0',
        ]
        tsv_lines = [gene_row('CM094067.1', 250000, 350000)]
        rows = self._run(tsv_lines, bed_lines, tsv_name='HG00097_hap1_chr13_intersect.tsv')
        self.assertEqual(rows[0][-1], 'TRUE')   # on_acrocentric_short_arm

    def test_acrocentric_arm_false_when_right_of_hor(self):
        bed_lines = [
            'CM094067.1\t0\t200000\tactive_hor_1\t100\t.\t0\t200000\t255,0,0',
            'CM094067.1\t400000\t700000\tHSAT3\t100\t.\t400000\t700000\t0,0,255',
        ]
        tsv_lines = [gene_row('CM094067.1', 250000, 350000)]
        rows = self._run(tsv_lines, bed_lines, tsv_name='HG00097_hap1_chr13_intersect.tsv')
        self.assertEqual(rows[0][-1], 'FALSE')

    def test_non_acrocentric_chrom_always_false_for_arm(self):
        bed_lines = [
            'CM094067.1\t0\t200000\tactive_hor_1\t100\t.\t0\t200000\t255,0,0',
        ]
        tsv_lines = [gene_row('CM094067.1', 50000, 80000)]
        rows = self._run(tsv_lines, bed_lines, tsv_name='HG00097_hap1_chr1_intersect.tsv')
        self.assertEqual(rows[0][-1], 'FALSE')

    def test_between_true_after_merging_small_adjacent_blocks(self):
        # two small non-CT blocks (each 60 kb, individually < 100 kb) that are
        # adjacent — after merging they become 120 kb → should pass the threshold
        bed_lines = [
            'CM094067.1\t0\t60000\tHSAT2\t100\t.\t0\t60000\t0,0,255',
            'CM094067.1\t60000\t120000\tHSAT3\t100\t.\t60000\t120000\t0,0,255',
            'CM094067.1\t120000\t300000\tct_1_1\t100\t.\t120000\t300000\t224,224,224',
            'CM094067.1\t300000\t360000\tmon\t100\t.\t300000\t360000\t0,255,0',
            'CM094067.1\t360000\t500000\tdhor\t100\t.\t360000\t500000\t0,255,0',
        ]
        tsv_lines = [gene_row('CM094067.1', 150000, 250000)]
        rows = self._run(tsv_lines, bed_lines)
        self.assertEqual(rows[0][-2], 'TRUE')  # merged left=120kb, merged right=200kb

    def test_between_true_when_non_nearest_large_blocks_present(self):
        # nearest neighbors are small (50kb each, < 100kb), but farther blocks are > 100kb
        # → new logic should return TRUE
        bed_lines = [
            'CM094067.1\t0\t150000\tHSAT2\t100\t.\t0\t150000\t0,0,255',           # far left: 150kb
            'CM094067.1\t150000\t200000\tHSAT3\t100\t.\t150000\t200000\t0,0,255',  # near left: 50kb
            'CM094067.1\t200000\t400000\tct_1_1\t100\t.\t200000\t400000\t224,224,224',
            'CM094067.1\t400000\t450000\tmon\t100\t.\t400000\t450000\t0,255,0',    # near right: 50kb
            'CM094067.1\t450000\t700000\tdhor\t100\t.\t450000\t700000\t0,255,0',   # far right: 250kb
        ]
        tsv_lines = [gene_row('CM094067.1', 250000, 350000)]
        rows = self._run(tsv_lines, bed_lines)
        self.assertEqual(rows[0][-2], 'TRUE')

    def test_short_row_gets_false_false(self):
        bed_lines = ['chr1\t0\t200000\tHSAT2\t100\t.\t0\t200000\t0,0,255']
        tsv_lines = ['only_two_cols\tsomething']
        rows = self._run(tsv_lines, bed_lines)
        self.assertEqual(rows[0][-2], 'FALSE')
        self.assertEqual(rows[0][-1], 'FALSE')

    def test_all_acrocentric_chroms_recognized(self):
        for chrom in ['chr13', 'chr14', 'chr15', 'chr21', 'chr22']:
            self.assertIn(chrom, ACROCENTRIC)

    def test_non_acrocentric_not_in_set(self):
        for chrom in ['chr1', 'chr2', 'chrX', 'chrY']:
            self.assertNotIn(chrom, ACROCENTRIC)


if __name__ == '__main__':
    unittest.main(verbosity=2)
