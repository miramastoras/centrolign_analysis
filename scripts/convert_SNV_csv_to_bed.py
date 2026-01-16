#!/usr/bin/env python3

import pandas as pd
import glob
import os
import argparse


# ------------------ ASAT HELPERS ------------------

def load_bed_files(bed_folder, chrom):
    """
    Load *_asat_arrays.bed files and return:
        { sample_name: { "contig": contig, "start": start } }
    """
    bed_dict = {}

    for fname in os.listdir(bed_folder):
        if not fname.endswith("_asat_arrays.bed"):
            continue

        sample = fname.replace("_asat_arrays.bed", "")
        path = os.path.join(bed_folder, fname)

        with open(path) as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue

                contig, start, end, chrom_field = line.rstrip("\n").split("\t")[:4]
                if chrom_field == chrom:
                    bed_dict[sample] = {
                        "contig": contig,
                        "start": int(start)
                    }
                    break

    return bed_dict


# ------------------ MAIN ------------------

def main():
    parser = argparse.ArgumentParser(
        description="Convert SNV CSVs into a single BED with query and reference rows."
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Folder with SNV CSV files")
    parser.add_argument("-o", "--output", required=True,
                        help="Output folder for BED files")
    parser.add_argument("-a", "--asat_beds", required=True,
                        help="Folder with *_asat_arrays.bed files")
    parser.add_argument("-c", "--chrom", required=True,
                        help="Chromosome to select ASAT regions (e.g. chr1)")

    args = parser.parse_args()

    input_dir = os.path.abspath(args.input)
    output_dir = os.path.abspath(args.output)
    asat_dir = os.path.abspath(args.asat_beds)

    os.makedirs(output_dir, exist_ok=True)

    asat = load_bed_files(asat_dir, args.chrom)
    if not asat:
        raise RuntimeError("No ASAT entries found for requested chromosome")

    csv_files = glob.glob(os.path.join(input_dir, "*.csv"))
    if not csv_files:
        raise FileNotFoundError("No CSV files found")

    for csv_path in csv_files:
        df = pd.read_csv(csv_path)

        ref_id = df["ref_id"].iloc[0]
        qry_id = df["qry_id"].iloc[0]

        if ref_id not in asat:
            raise KeyError(f"Reference '{ref_id}' missing from ASAT beds")
        if qry_id not in asat:
            raise KeyError(f"Query '{qry_id}' missing from ASAT beds")

        ref_contig = asat[ref_id]["contig"]
        ref_offset = asat[ref_id]["start"]

        qry_contig = asat[qry_id]["contig"]
        qry_offset = asat[qry_id]["start"]

        # Reference rows
        ref_df = pd.DataFrame({
            "contig": ref_contig,
            "start": df["ref_pos"] + ref_offset,
            "end": df["ref_pos"] + ref_offset + 1 ,
            "role": "reference",
            "local_snv_rate_percentile": df["local_snv_rate_percentile"],
            "dist_to_break": df["dist_to_break"],
        })

        # Query rows
        qry_df = pd.DataFrame({
            "contig": qry_contig,
            "start": df["qry_pos"] + qry_offset,
            "end": df["qry_pos"] + qry_offset +1,
            "role": "query",
            "local_snv_rate_percentile": df["local_snv_rate_percentile"],
            "dist_to_break": df["dist_to_break"],
        })

        out_df = pd.concat([ref_df, qry_df], ignore_index=True)

        # Get input filename (no directories)
        csv_name = os.path.basename(csv_path)

        # Replace .csv with .bed
        bed_name = os.path.splitext(csv_name)[0] + ".bed"

        out_path = os.path.join(output_dir, bed_name)

        out_df.to_csv(out_path, sep="\t", header=False, index=False)


if __name__ == "__main__":
    main()
