## Deriving SV calls from induced pairwise cigar strings

### Step 1: Identify low divergence clades to call variants in

```sh
python3 /private/home/mmastora/github/centrolign_analysis/scripts/identify_low_divergence_clades.py \
  --newick /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr8/split_by_2/subgroup_1_tree.nwk \
  --distance_csv /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/chr8_r2_QC_v2_centrolign_pairwise_distance.csv \
  --output_prefix /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/subgroup_1
```


### Step 2: Write python script to call SVs from cigar strings

input: directory containing cigar strings
output: bedPE file with SV coords relative to each assembly

```sh
python3  /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/call_SVs_pairwise.py -c /Users/miramastoras/Desktop/chr12_cigars/pairwise_cigar_ > /Users/miramastoras/Desktop/SV_bedfile.bed
```

```
def plot_SV_length_dist(sv_positions):
    """
    Plot a swarm plot showing the length distribution of structural variants (SVs)
    by type (Insertion or Deletion).

    Parameters
    ----------
    sv_positions : list of tuples
        Each tuple should have the format (ref, start_ref, end_ref, query, start_query, end_query, type),
        where type is either 'I' for insertion or 'D' for deletion.
    """
    # Compute SV lengths
    sv_data = []
    for ref, start_ref, end_ref, query, start_query, end_query, sv_type in sv_positions:
        if sv_type == 'I':
            length = abs(end_query - start_query)
        elif sv_type == 'D':
            length = abs(end_ref - start_ref)
        else:
            continue  # skip unknown types

        sv_data.append({"Type": "Insertion" if sv_type == "I" else "Deletion",
                        "Length": length})

    # Convert to DataFrame for easy plotting
    df = pd.DataFrame(sv_data)
    print(df)
    # Plot
    plt.figure(figsize=(8, 6))
    sns.violinplot(x="Type", y="Length", data=df, size=6, palette="Set2")
    plt.title("Structural Variant Length Distribution", fontsize=14)
    plt.ylabel("SV Length (bp)")
    plt.xlabel("SV Type")
    plt.tight_layout()
    plt.savefig('/Users/miramastoras/Desktop/svs.png')
```
