#!/usr/bin/env python3
"""
Unit tests for maf_to_bed_per_sample.py

Run with:
    python test_maf_to_bed_per_sample.py
"""

import os
import sys
import tempfile
import traceback

from maf_to_bed_per_sample import sample_name, chrom_name, process_block


# ── Helpers ───────────────────────────────────────────────────────────────────

def make_state(outdir):
    """Return (seen_samples, write_record) closure as used in main()."""
    seen = set()
    def write_record(samp, record):
        mode = "a" if samp in seen else "w"
        seen.add(samp)
        with open(os.path.join(outdir, f"{samp}.bed"), mode) as fh:
            fh.write(record)
    return seen, write_record


def read_bed(path):
    rows = []
    with open(path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            rows.append((fields[0], int(fields[1]), int(fields[2]), fields[3], int(fields[4])))
    return rows


# ── Test runner ───────────────────────────────────────────────────────────────

TESTS = []

def test(fn):
    TESTS.append(fn)
    return fn

def run_all():
    passed = failed = 0
    for fn in TESTS:
        try:
            fn()
            print(f"  PASS  {fn.__name__}")
            passed += 1
        except Exception:
            print(f"  FAIL  {fn.__name__}")
            traceback.print_exc()
            failed += 1
    print(f"\n{passed}/{passed+failed} tests passed")
    return failed


# ── sample_name ───────────────────────────────────────────────────────────────

@test
def test_sample_name_two_parts():
    assert sample_name("HG00097.1") == "HG00097.1"

@test
def test_sample_name_many_parts():
    assert sample_name("HG00097.1.CM094066.1") == "HG00097.1"

@test
def test_sample_name_hap2():
    assert sample_name("HG00097.2.CM094066.1") == "HG00097.2"

@test
def test_sample_name_no_dot():
    assert sample_name("CHM13") == "CHM13"


# ── chrom_name ────────────────────────────────────────────────────────────────

@test
def test_chrom_name_standard():
    assert chrom_name("CHM13.chr10") == "chr10"

@test
def test_chrom_name_no_dot():
    assert chrom_name("chr10") == "chr10"


# ── process_block ─────────────────────────────────────────────────────────────

@test
def test_basic_block():
    with tempfile.TemporaryDirectory() as outdir:
        _, write_record = make_state(outdir)
        block = [
            "s CHM13.chr10        100  10  +  134758134  ACGTACGTAC".split(),
            "s HG00097.1.chr10    200  10  +  134577915  ACGTACGTAC".split(),
        ]
        process_block(block, "CHM13", write_record)
        rows = read_bed(os.path.join(outdir, "HG00097.1.bed"))
        assert len(rows) == 1
        chrom, start, end, samp, gaps = rows[0]
        assert chrom == "chr10"
        assert start == 100
        assert end   == 110
        assert samp  == "HG00097.1"
        assert gaps  == 0

@test
def test_gap_counting():
    with tempfile.TemporaryDirectory() as outdir:
        _, write_record = make_state(outdir)
        block = [
            "s CHM13.chr10     100  10  +  134758134  ACGTACGTAC".split(),
            "s HG00097.1.chr10 200   7  +  134577915  ACG---GTAC".split(),
        ]
        process_block(block, "CHM13", write_record)
        rows = read_bed(os.path.join(outdir, "HG00097.1.bed"))
        assert rows[0][4] == 3

@test
def test_insertion_columns_excluded():
    with tempfile.TemporaryDirectory() as outdir:
        _, write_record = make_state(outdir)
        # '--' in reference = insertion in sample; should not affect ref coords or gap count
        block = [
            "s CHM13.chr10     100  10  +  134758134  ACGT--ACGTAC".split(),
            "s HG00097.1.chr10 200  12  +  134577915  ACGTTTACGTAC".split(),
        ]
        process_block(block, "CHM13", write_record)
        rows = read_bed(os.path.join(outdir, "HG00097.1.bed"))
        assert rows[0][1] == 100
        assert rows[0][2] == 110
        assert rows[0][4] == 0

@test
def test_two_samples():
    with tempfile.TemporaryDirectory() as outdir:
        _, write_record = make_state(outdir)
        block = [
            "s CHM13.chr10        100  10  +  134758134  ACGTACGTAC".split(),
            "s HG00097.1.chr10    200  10  +  100000000  ACGTACGTAC".split(),
            "s HG00438.2.chr10    300  10  +  100000000  ACGT--GTAC".split(),
        ]
        process_block(block, "CHM13", write_record)
        assert os.path.exists(os.path.join(outdir, "HG00097.1.bed"))
        assert os.path.exists(os.path.join(outdir, "HG00438.2.bed"))
        assert read_bed(os.path.join(outdir, "HG00097.1.bed"))[0][4] == 0
        assert read_bed(os.path.join(outdir, "HG00438.2.bed"))[0][4] == 2

@test
def test_no_ref_block_skipped():
    with tempfile.TemporaryDirectory() as outdir:
        _, write_record = make_state(outdir)
        block = [
            "s HG00097.1.chr10  200  10  +  100000000  ACGTACGTAC".split(),
            "s HG00438.2.chr10  300  10  +  100000000  ACGTACGTAC".split(),
        ]
        process_block(block, "CHM13", write_record)
        assert len(os.listdir(outdir)) == 0

@test
def test_multiple_blocks_appended():
    with tempfile.TemporaryDirectory() as outdir:
        seen, write_record = make_state(outdir)
        block1 = [
            "s CHM13.chr10        100  10  +  134758134  ACGTACGTAC".split(),
            "s HG00097.1.chr10    200  10  +  100000000  ACGTACGTAC".split(),
        ]
        block2 = [
            "s CHM13.chr10        500   8  +  134758134  GGGGCCCC".split(),
            "s HG00097.1.chr10    600   8  +  100000000  GGGG--CC".split(),
        ]
        process_block(block1, "CHM13", write_record)
        process_block(block2, "CHM13", write_record)
        rows = read_bed(os.path.join(outdir, "HG00097.1.bed"))
        assert len(rows) == 2
        assert rows[0] == ("chr10", 100, 110, "HG00097.1", 0)
        assert rows[1] == ("chr10", 500, 508, "HG00097.1", 2)

@test
def test_minigraph_skipped():
    with tempfile.TemporaryDirectory() as outdir:
        _, write_record = make_state(outdir)
        block = [
            "s CHM13.chr10          100  10  +  134758134  ACGTACGTAC".split(),
            "s _MINIGRAPH_.s48827     0  10  +      10000  ACGTACGTAC".split(),
            "s HG00097.1.chr10      200  10  +  100000000  ACGTACGTAC".split(),
        ]
        process_block(block, "CHM13", write_record)
        assert not os.path.exists(os.path.join(outdir, "_MINIGRAPH_.bed"))
        assert os.path.exists(os.path.join(outdir, "HG00097.1.bed"))

@test
def test_minus_strand_ref_coords():
    with tempfile.TemporaryDirectory() as outdir:
        _, write_record = make_state(outdir)
        block = [
            "s CHM13.chr10        100  10  +  134758134  ACGTACGTAC".split(),
            "s HG00097.1.chr10    200  10  -  100000000  ACGTACGTAC".split(),
        ]
        process_block(block, "CHM13", write_record)
        rows = read_bed(os.path.join(outdir, "HG00097.1.bed"))
        assert rows[0][0] == "chr10"
        assert rows[0][1] == 100
        assert rows[0][2] == 110


if __name__ == "__main__":
    sys.exit(run_all())