#!/usr/bin/env python3
"""
Convert SV bed pe files to assembly coordinates.
"""

import argparse
import os

def process_bed_file(input_path, asat_dict, output_path, fmt):
    """
    Convert coordinates in BED or BEDPE files while preserving all columns.
    Replace contig names with contigs from asat BEDs.
    """

    with open(input_path, "r") as infile, open(output_path, "w") as outfile:
        for raw_line in infile:
            if raw_line.startswith("#") or not raw_line.strip():
                outfile.write(raw_line)
                continue

            cols = raw_line.rstrip("\n").split("\t")

            if fmt == "bed":
                if len(cols) < 3:
                    raise ValueError(f"Invalid BED line (<3 columns): {raw_line}")

                smp = cols[0]

                asat_contig, asat_start = asat_dict[smp][0], int(asat_dict[smp][1])

                # replace contig
                cols[0] = asat_contig

                # shift coordinates
                cols[1] = str(int(cols[1]) + asat_start)
                cols[2] = str(int(cols[2]) + asat_start)

            elif fmt == "bedpe":
                if len(cols) < 6:
                    raise ValueError(f"Invalid BEDPE line (<6 columns): {raw_line}")

                smp1 = cols[0]
                smp2 = cols[3]

                contig1, start1 = asat_dict[smp1][0], int(asat_dict[smp1][1])
                contig2, start2 = asat_dict[smp2][0], int(asat_dict[smp2][1])

                # replace contigs
                cols[0] = contig1
                cols[3] = contig2

                # shift coordinates
                cols[1] = str(int(cols[1]) + start1)
                cols[2] = str(int(cols[2]) + start1)
                cols[4] = str(int(cols[4]) + start2)
                cols[5] = str(int(cols[5]) + start2)

            outfile.write("\t".join(cols) + "\n")


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
        "--format",
        choices=["bed", "bedpe"],
        required=True,
        help="Input file format: bed or bedpe"
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
            output_path = os.path.join(output_folder, filename)

            process_bed_file(input_path, asat_dict, output_path, args.format)


if __name__ == "__main__":
    main()
