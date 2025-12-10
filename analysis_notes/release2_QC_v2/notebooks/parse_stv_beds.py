import glob
import re
import os
import numpy as np
import pandas as pd

def assign_chr_from_asat(input_df,
                         folder="/private/groups/migalab/juklucas/censat_regions/active_arrays/"):
    """
    Given an input dataframe with columns ["smp", "start", "end", "name"],
    read all ASAT array files from the folder and map smp → chrom_assignment.
    
    Returns the input dataframe with an added 'chr' column.
    """

    # -----------------------------
    # 1. Find all ASAT array TSV files
    # -----------------------------
    pattern = os.path.join(folder, "asat_arrays_*_raw.tsv")
    files = glob.glob(pattern)

    if not files:
        raise FileNotFoundError(f"No ASAT array files found in {folder}")

    lookup_tables = []

    # -----------------------------
    # 2. Read each file
    # -----------------------------
    for f in files:
        # Extract chromosome from filename
        # asat_arrays_chr1_raw.tsv → chr1
        df = pd.read_csv(f, sep="\t")

        # Keep only the needed columns
        if not {"sequence_id", "chrom_assignment", "assembly_id"} <= set(df.columns):
            raise ValueError(f"File {f} missing required columns.")

        lookup_tables.append(df[["sequence_id", "chrom_assignment", "assembly_id"]])

    # -----------------------------
    # 3. Combine all lookup tables
    # -----------------------------
    lookup_df = pd.concat(lookup_tables, ignore_index=True).drop_duplicates()

    # -----------------------------
    # 4. Merge to input dataframe
    # -----------------------------
    merged = input_df.merge(
        lookup_df,
        left_on="smp",
        right_on="sequence_id",
        how="left"
    )

    # Rename chrom_assignment to chr
    merged = merged.rename(columns={"chrom_assignment": "chr"})

    # Drop redundant sequence_id column
    merged = merged.drop(columns=["sequence_id"])

    return merged


def aggregate_hor_stats(pattern):
    files = glob.glob(pattern)
    all_dfs = []

    for f in files:
        # Extract assembly_id from path
        # Example path: .../STV/<assembly_id>/analysis/...
        parts = f.split("/")
        assembly_id = parts[parts.index("STV") + 1]

        # Load file, skip first 2 header lines
        df = pd.read_csv(f, skiprows=2, names=["name", "count"], sep="\t")
        # Add assembly_id as a new column
        df["assembly_id"] = assembly_id

        all_dfs.append(df)

    # Combine everything
    big_df = pd.concat(all_dfs, ignore_index=True)

    return big_df


def load_beds(pattern_or_path):
    """
    Accepts either:
      - a single BED file path
      - a glob pattern that matches many BED files

    Returns a concatenated dataframe with columns:
        start, end, name, length
    """
    # Check if it's a single file
    if os.path.isfile(pattern_or_path):
        bed_files = [pattern_or_path]
    else:
        bed_files = glob.glob(pattern_or_path, recursive=True)

    if len(bed_files) == 0:
        raise FileNotFoundError(f"No BED files matched: {pattern_or_path}")

    dfs = []
    for bed in bed_files:
        df = pd.read_csv(
            bed,
            sep="\t",
            header=None,
            usecols=[0,1, 2, 3],  # start, end, name
            names=["smp","start", "end", "name"]
        )
        df["length"] = df["end"] - df["start"]
        dfs.append(df)

    return pd.concat(dfs, ignore_index=True)


def print_name_range_stats(df, freq_df):
    """
    For each 'name' with frequency > 5%, print:
      - the sample (smp) with the maximum range (end - start)
      - the sample (smp) with the minimum range (end - start)
    
    Parameters
    ----------
    df : DataFrame
        Original BED-like dataframe with columns start, end, smp, name
    freq_df : DataFrame
        DataFrame with columns 'name' and 'freq', filtered for freq > 0.05
    """

    # Ensure numeric coordinates
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"]   = pd.to_numeric(df["end"], errors="coerce")

    # Compute range
    df["range"] = df["end"] - df["start"]

    # Only keep names with freq > 5%
    valid_names = set(freq_df["name"])
    df = df[df["name"].isin(valid_names)]

    for name, group in df.groupby("name"):
        max_row = group.loc[group["range"].idxmax()]
        min_row = group.loc[group["range"].idxmin()]

        print(f"\n{name}:")
        print(f"  Max range = {max_row['range']} bp "
              f"(start={max_row['start']}, end={max_row['end']}, sample={max_row['smp']})")
        print(f"  Min range = {min_row['range']} bp "
              f"(start={min_row['start']}, end={min_row['end']}, sample={min_row['smp']})")

def summarize_beds_per_chr(df):
    """
    Given a dataframe with columns:
    smp, start, end, name, length, chr, assembly_id

    For each chromosome + name:
        - compute min_range and max_range (max can only come from rows with range <= 2*min_range) to avoid rows where the adjacent STVs are merged.
        - compute a frequency score:
              if range <= 2*min_range → +1
              else → +(range / mean(min_range, max_range)). For merged windows we need to divide by the average length of the HOR to 
              get the correct count of HOR STVs in that window
    Returns a DataFrame with columns:
        chr, name, min_range, max_range, frequency
    """
    # Ensure numeric
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"]   = pd.to_numeric(df["end"], errors="coerce")

    # Compute range
    df["range"] = df["end"] - df["start"]

    summary = []

    # Group by chromosome + name
    for (chrom, name), group in df.groupby(["chr", "name"]):

        if group.empty:
            continue

        # min range
        min_range = group["range"].min()

        # rows allowed to define max_range
        filtered = group[group["range"] <= 2 * min_range]

        if filtered.empty:
            max_range = np.nan
        else:
            max_range = filtered["range"].max()

        # mean range for scaling
        mean_range = np.nanmean([min_range, max_range])

        # compute weighted counts
        counts = 0
        for r in group["range"]:
            if r <= 2 * min_range:
                counts += 1
            else:
                counts += r / mean_range

        summary.append({
            "chr": chrom,
            "name": name,
            "min_range": min_range,
            "max_range": max_range,
            "counts": counts
        })

    # Convert summary list to DataFrame
    summary_df = pd.DataFrame(summary)

    # --- Add per-chromosome frequency ---
    # total counts per chromosome
    chrom_totals = summary_df.groupby("chr")["counts"].transform("sum")

    summary_df["chrom_freq"] = summary_df["counts"] / chrom_totals

    return summary_df


def summarize_beds(df):
    """
    Given a dataframe with columns start, end, name, length,
    compute counts and the min start and max end per name.
    For max end, drop any rows where end > 2 * min start.
    """

    # Ensure numeric start/end
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"]   = pd.to_numeric(df["end"], errors="coerce")

    # Compute length/range
    df["range"] = df["end"] - df["start"]

    summary_list = []

    for name, group in df.groupby("name"):
        if group.empty:
            continue

        # Compute min range
        min_range_val = group["range"].min()

        # Filter rows for max range
        filtered_group = group[group["range"] <= 2 * min_range_val]

        # Compute max range
        max_range_val = filtered_group["range"].max() if not filtered_group.empty else float("nan")

        # Count total rows
        count_val = len(group)

        summary_list.append({
            "name": name,
            "count": count_val,
            "min_range": min_range_val,
            "max_range": max_range_val
        })

    return pd.DataFrame(summary_list)

def extract_chr(name):
    """
    Extract chromosome from 'name' column safely.

    Matches:
        C1–C22
        CX, CY

    Returns 'chr#'
    """
    m = re.search(r"C(\d{1,2}|X|Y)", name)
    if m:
        return f"chr{m.group(1)}"
    return None


def compute_frequencies(df):
    """
    Adds:
        - chrom
        - freq
    and filters to freq > 5%
    """
    df["chrom"] = df["name"].apply(extract_chr)

    chrom_totals = df.groupby("chrom")["count"].transform("sum")
    df["freq"] = df["count"] / chrom_totals

    return df[df["freq"] > 0.05].copy()


# -------------------------------
# Full pipeline call
# -------------------------------

pattern = "/private/groups/migalab/hloucks/DeepLineage/DeepLineage_censat/STV/*/analysis/HOR-StV_outputs/*/*HOR_StV_row.bed"
#pattern="/private/groups/patenlab/mira/centrolign/annotations/CHM13v2.0_stv_raw.bed"

# Load + merge one or many files
big_df = load_beds(pattern)

# assign chromosomes based on julian's QC mapping 
big_df_chroms=assign_chr_from_asat(big_df)

results_df=summarize_beds_per_chr(big_df_chroms)

filtered = results_df[results_df["chrom_freq"] > 0.05]

# write to file 
filtered.to_csv("/private/groups/patenlab/mira/centrolign/annotations/stvs_gt.05_12092025_chrfix.tsv", sep="\t", index=False)
#filter to frequencies > 0.05 


# # Compute counts, mean ranges
#result_df = summarize_beds(big_df)

# # Add chromosome + frequency, filter to freq > 5%
# filtered_df = compute_frequencies(result_df)

# # Print min/max ranges only for names with freq > 5%
# #print_name_range_stats(big_df, filtered_df)

# # Write it out
# outfile = "/private/groups/patenlab/mira/centrolign/annotations/stvs_gt.05_12082025.tsv"
# #outfile = "/private/groups/patenlab/mira/centrolign/annotations/CHM13_stvs_gt.05.tsv"
# filtered_df.to_csv(outfile, sep="\t", index=False)

# print("Wrote:", outfile)
# #print(filtered_df.head())

# ## also write out all frequencies to check 

# big_df["chrom"] = big_df["name"].apply(extract_chr)

# # List of chromosomes you want to keep
# chroms_to_keep = ["chr22", "chr19", "chr5", "21"]

# # Subset the DataFrame
# subset_df = big_df[big_df["chrom"].isin(chroms_to_keep)]

# print(subset_df.head())
# #big_df.to_csv("/private/groups/patenlab/mira/centrolign/annotations/all_stvs_12092025.tsv", sep="\t", index=False)
