## Deriving SV calls from induced pairwise cigar strings

### Starting with Chr8 as a test case

#### Step 1: Identify low divergence clades to call variants in

```sh
docker run -u `id -u`:`id -g` -v /private/groups:/private/groups/ \
    miramastoras/centromere_scripts:v0.1.4 \
    python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/identify_low_divergence_clades.py \
  --newick /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr8/split_by_2/subgroup_1_tree.nwk \
  --distance_csv /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/chr8_r2_QC_v2_centrolign_pairwise_distance.csv \
  --max_pairwise_dist 0.8 \
  --output_prefix /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/subgroup_1
```
chr8 subgroup 1 - all one clade

chr 8 subgroup 0 - 8 clades
```sh
docker run -u `id -u`:`id -g` -v /private/groups:/private/groups/ \
    miramastoras/centromere_scripts:v0.1.4 \
    python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/identify_low_divergence_clades.py \
  --newick /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr8/split_by_2/subgroup_0_tree.nwk \
  --distance_csv /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/chr8_r2_QC_v2_centrolign_pairwise_distance.csv \
  --max_pairwise_dist 0.8 \
  --output_prefix /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/subgroup_0
```

```sh
chromosomes=("chr8")

for chr in "${chromosomes[@]}"
do
  python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap_v2.py \
    -t /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_all_pairs_nj_tree.format5.nwk \
    -s /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}.samples.txt \
    -p /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_pairwise_distance.csv \
    -m "Centrolign all pairs distances" \
    -n "${chr} NJ tree" \
    -d "All pairs Distances" \
    -o /Users/miramastoras/Desktop/color_subgroups_heatmap/${chr}_r2_QC_v2_clade_2_subgroup0 --no_labels \
    --highlight_samples /Users/miramastoras/Desktop/color_subgroups_heatmap/chr8_subgroup0_clade2.txt
done
```

#### Step 2: Call SVs from cigar strings within each clade

input: directory containing cigar strings
output: bedPE file with SV coords relative to each assembly

Testing how long it takes to run on chr8 subgroup 1
```sh
time python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/call_SVs_pairwise.py \
  -c /Users/miramastoras/Desktop/chr12_cigars/pairwise_cigar_ \
  -s /Users/miramastoras/Desktop/chr12_cigars/samples.txt \
  -o /Users/miramastoras/Desktop/chr12_cigars/SVs_
```

```sh
time python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/call_SVs_pairwise.py \
  -c /Users/miramastoras/Desktop/chr12_cigars/pairwise_cigar_ \
  -s /Users/miramastoras/Desktop/chr12_cigars/samples.txt \
  -o /Users/miramastoras/Desktop/chr12_cigars/SVs_
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
