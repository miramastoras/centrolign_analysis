#!/usr/bin/env python3
"""
Convert SV bed pe files to assembly coordinates. 
"""

import argparse
import os

def process_bedpe_file(input_path, asat_dict, output_path):
    """Reads a BEDPE file line by line and writes it to a new file."""
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            line = line.strip().split("\t")
            smp1 = line[0]
            smp2 = line[3]

            # start coord of asat array for each sample
            smp1_offset = int(asat_dict[smp1][1])
            smp2_offset = int(asat_dict[smp2][1])

            # print contig instead of asm name
            smp1_contig = asat_dict[smp1][0]
            smp2_contig = asat_dict[smp2][0]

            # add asat start coord to convert to asm coords
            smp1_start = int(line[1])+ smp1_offset
            smp1_end = int(line[2]) + smp1_offset

            smp2_start = int(line[4]) + smp2_offset
            smp2_end = int(line[5]) + smp2_offset

            print(smp1_contig, smp1_start, smp1_end, smp2_contig, smp2_start, smp2_end, line[6], line[7], sep="\t", file=outfile)

def load_bed_files(bed_folder, chrom):
    """
    Reads BED files from bed_folder, filters lines where chrom matches,
    and stores them in a dictionary:
        { sample_name: [contig, start, end, chrom] }

    If multiple lines match the chromosome, only the first match is stored.
    """
    bed_dict = {}

    for filename in os.listdir(bed_folder):
        if not filename.endswith("_asat_arrays.bed"):
            continue

        sample_name = filename.replace("_asat_arrays.bed", "")
        file_path = os.path.join(bed_folder, filename)

        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                cols = line.split('\t')
                if len(cols) < 4:
                    continue  # skip malformed lines

                contig, start, end, chrom_field = cols[:4]
                if chrom_field == chrom:
                    bed_dict[sample_name] = [contig, start, end, chrom_field]
                    break  # stop after first match for this sample

    return bed_dict


def main():
    parser = argparse.ArgumentParser(
        description="Process BEDPE and BED files, filtering BED entries by chromosome."
    )
    parser.add_argument(
        "-s", "--sv_beds",
        required=True,
        help="Folder containing BEDPE files (sample1_sample2.bed)."
    )
    parser.add_argument(
        "-a", "--asat_beds",
        required=True,
        help="Folder containing asat annotation BED files (sample_asat_arrays.bed)."
    )
    parser.add_argument(
        "-c", "--chrom",
        required=True,
        help="Chromosome to filter BED files by (e.g., chr12)."
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output folder for processed BEDPE files."
    )

    args = parser.parse_args()

    sv_beds_folder = os.path.abspath(args.sv_beds)
    asat_beds_folder = os.path.abspath(args.asat_beds)
    output_folder = os.path.abspath(args.output)

    # Validate folders
    if not os.path.isdir(sv_beds_folder):
        raise NotADirectoryError(f"BEDPE folder not found: {sv_beds_folder}")
    if not os.path.isdir(asat_beds_folder):
        raise NotADirectoryError(f"BED folder not found: {asat_beds_folder}")

    os.makedirs(output_folder, exist_ok=True)

    # Load Asat bed files into a dict
    asat_dict = load_bed_files(asat_beds_folder, args.chrom)

    # Process BEDPE files
    for filename in os.listdir(sv_beds_folder):
        if filename.endswith(".bed"):
            input_path = os.path.join(sv_beds_folder, filename)
            output_filename = filename.replace(".bed", "_asm_coords.bed")
            output_path = os.path.join(output_folder, output_filename)

            process_bedpe_file(input_path, asat_dict,output_path)

if __name__ == "__main__":
    main()
