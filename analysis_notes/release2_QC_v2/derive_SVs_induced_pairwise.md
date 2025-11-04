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
cut -f1 subgroup_1_seqs.fasta.fai  > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/chr8_subgroup1_samples.txt
```

```sh
#!/bin/bash
#SBATCH --job-name=SVs_chr8_subgroup_1
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=call_SVs_%x.%j.log
#SBATCH --time=7-00:00

time python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/call_SVs_pairwise.py \
  -c /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr8/subgroup_1/induced_pairwise_cigars/pairwise_cigar_ \
  -s /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/chr8_subgroup1_samples.txt \
  -o /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/subgroup_1_SV_beds/
```
Ran in two minutes

Local computer testing:
```sh
time python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/call_SVs_pairwise.py \
  -c /Users/miramastoras/Desktop/chr12_cigars/pairwise_cigar_ \
  -s /Users/miramastoras/Desktop/chr12_cigars/samples.txt \
  -o /Users/miramastoras/Desktop/chr12_cigars/SVs_
```
Plot SV length distributions
```sh
docker run -u `id -u`:`id -g` -v /private/groups:/private/groups/ \
    miramastoras/centromere_scripts:v0.1.4 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/plot_SVs_pairwise.py \
    -b /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/subgroup_1_SV_beds/ \
    -o /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/SV_violin_subgroup1.png
```
