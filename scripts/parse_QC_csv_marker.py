import argparse
import pandas as pd
import os

def main():
    parser = argparse.ArgumentParser(
        description="Parse HGSVC asat arrays CSV and output per-sample BED files")
    parser.add_argument("qc_table", help="Path to hgsvc_asat_arrays.csv")
    parser.add_argument("output_path", help="Output directory for BED files")
    args = parser.parse_args()

    df = pd.read_csv(args.qc_table)

    os.makedirs(args.output_path, exist_ok=True)

    for assembly_id, group in df.groupby("assembly_id"):
        sample_id = group["assembly"].iloc[0]
        haplotype  = group["haplotype"].iloc[0]

        # For each (sequence_id, chrom_assignment), merge all asat intervals
        # into a single spanning region: min(asat_start) to max(asat_end)
        merged = (
            group.groupby(["sequence_id", "chrom_assignment"], as_index=False)
            .agg(asat_start=("asat_start", "min"),
                 asat_end=("asat_end",   "max"))
            [["sequence_id", "asat_start", "asat_end", "chrom_assignment"]]
        )

        filename = os.path.join(args.output_path, f"{sample_id}.{haplotype}_asat_arrays.bed")
        file_exists = os.path.exists(filename)

        merged.to_csv(
            filename,
            sep="\t",
            header=False,
            index=False,
            mode='a' if file_exists else 'w'
        )

if __name__ == "__main__":
    main()
