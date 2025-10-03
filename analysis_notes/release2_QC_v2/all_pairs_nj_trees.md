## Create neighbor-joining trees based on all pairs alignments, and decide where to split trees for MSA runs

### Chr 1, Chr 12, Chr 10, Chr 11, Chr 8

Get distance matrices for all pairs
```sh
chromosomes=("chr1" "chr12" "chr10" "chr11" "chr8")

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices

git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/ pull

for chr in "${chromosomes[@]}"
do
    python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/cigar_to_distance.py -d 1 \
    -a /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/${chr}/pairwise_cigar/ \
    -o /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/${chr}_r2_QC_v2_centrolign_
  done
```

Generate NJ trees for all samples (except chr6 and 12 since those have been created)
```sh
chromosomes=("chr1" "chr12" "chr10" "chr11" "chr8")

for chr in "${chromosomes[@]}"
do
    docker run -u `id -u`:`id -g` -v /private/groups:/private/groups/ \
        miramastoras/centromere_scripts:v0.1.4 \
        python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/data_exploration/infer_tree.py \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/${chr}/pairwise_cigar/ \
        > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/${chr}_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk

    docker run -u `id -u`:`id -g` -v /private/groups:/private/groups/ \
        miramastoras/centromere_scripts:v0.1.4 \
        python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/convert_tree_format.py \
        -t /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/${chr}_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk
  done
```
Create sample lists for the heatmaps
```sh
chromosomes=("chr1" "chr12" "chr10" "chr11" "chr8")

for chr in "${chromosomes[@]}"
do
  cut -f1-2 /private/groups/migalab/juklucas/censat_regions/active_arrays/asat_arrays_${chr}.tsv | grep -v "sample_id" | sed 's/\t/./g' > /private/groups/patenlab/mira/centrolign/analysis/HPRC_release2_QCv2/sample_lists/${chr}.samples.txt
done

ls | while read line ; do wc -l $line;done
```

Plot pairwise heatmaps
```sh
chromosomes=("chr1" "chr12" "chr10" "chr11" "chr8")

for chr in "${chromosomes[@]}"
do
    python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap_v2.py \
        -t /Users/miramastoras/Desktop/HPRCr2_QC_v2/${chr}_r2_QC_v2_centrolign_all_pairs_nj_tree.format5.nwk \
        -s /Users/miramastoras/Desktop/HPRCr2_QC_v2/${chr}.samples.txt \
        -p /Users/miramastoras/Desktop/HPRCr2_QC_v2/${chr}_r2_QC_v2_centrolign_pairwise_distance.csv \
        -m "Centrolign all pairs distances" \
        -n "${chr} NJ tree" \
        -d "All pairs Distances" \
        -o /Users/miramastoras/Desktop/github_repos/centrolign_analysis/analysis_notes/release2_QC_v2/figures/${chr}_r2_QC_v2_all_pairs --no_labels
done
```

### Chr 2, Chr 6, Chr 7, chr 5, chr9

Get distance matrices and create NJ trees
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/ pull

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/logs

sbatch \
    --job-name=chr2_r2_QCv2_tree \
    --export=CHR=chr2 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs_nj_tree.sh

sbatch \
    --job-name=chr6_r2_QCv2_tree \
    --export=CHR=chr6 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs_nj_tree.sh

sbatch \
    --job-name=chr7_r2_QCv2_tree \
    --export=CHR=chr7 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs_nj_tree.sh

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees

sbatch \
    --job-name=chr5_r2_QCv2_tree \
    --export=CHR=chr5 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs_nj_tree.sh


sbatch \
    --job-name=chr9_r2_QCv2_tree \
    --export=CHR=chr9 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs_nj_tree.sh
```

Create sample lists for the heatmaps
```sh
chromosomes=("chr2" "chr6" "chr7" "chr5" "chr9")

for chr in "${chromosomes[@]}"
do
  cut -f1-2 /private/groups/migalab/juklucas/censat_regions/active_arrays/asat_arrays_${chr}.tsv | grep -v "sample_id" | sed 's/\t/./g' > /private/groups/patenlab/mira/centrolign/analysis/HPRC_release2_QCv2/sample_lists/${chr}.samples.txt
done

# Sanity check the sample list is the correct size
cd /private/groups/patenlab/mira/centrolign/analysis/HPRC_release2_QCv2/sample_lists/
ls  | while read line ; do wc -l $line;done
```

Copy scripts for plots to local
```sh
mkdir -p /private/groups/patenlab/mira/HPRC_release2_QCv2_all_pairs_heatmaps/
cp /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/*format5* /private/groups/patenlab/mira/HPRC_release2_QCv2_all_pairs_heatmaps/

cp /private/groups/patenlab/mira/centrolign/analysis/HPRC_release2_QCv2/sample_lists/*samples.txt /private/groups/patenlab/mira/HPRC_release2_QCv2_all_pairs_heatmaps/

cp /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/* /private/groups/patenlab/mira/HPRC_release2_QCv2_all_pairs_heatmaps/
```

Plot pairwise heatmaps
```sh
chromosomes=("chr2" "chr6" "chr7" "chr5" "chr9")

for chr in "${chromosomes[@]}"
do
    python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap_v2.py \
        -t /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_all_pairs_nj_tree.format5.nwk \
        -s /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}.samples.txt \
        -p /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_pairwise_distance.csv \
        -m "Centrolign all pairs distances" \
        -n "${chr} NJ tree" \
        -d "All pairs Distances" \
        -o /Users/miramastoras/Desktop/github_repos/centrolign_analysis/analysis_notes/release2_QC_v2/figures/${chr}_r2_QC_v2_all_pairs --no_labels
done
```
### Chr 13, Chr 14, Chr 15, chr 16, chr17, chr18

Get distance matrices and create NJ trees
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/ pull

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/logs

sbatch \
    --job-name=chr13_r2_QCv2_tree \
    --export=CHR=chr13 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs_nj_tree.sh

sbatch \
    --job-name=chr14_r2_QCv2_tree \
    --export=CHR=chr14 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs_nj_tree.sh


sbatch \
    --job-name=chr17_r2_QCv2_tree \
    --export=CHR=chr17 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs_nj_tree.sh

sbatch \
    --job-name=chr18_r2_QCv2_tree \
    --export=CHR=chr18 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs_nj_tree.sh
```
Create sample lists for the heatmaps
```sh
chromosomes=("chr13" "chr14" "chr17" "chr18")

for chr in "${chromosomes[@]}"
do
  cut -f1-2 /private/groups/migalab/juklucas/censat_regions/active_arrays/asat_arrays_${chr}.tsv | grep -v "sample_id" | sed 's/\t/./g' > /private/groups/patenlab/mira/centrolign/analysis/HPRC_release2_QCv2/sample_lists/${chr}.samples.txt
done

# Sanity check the sample list is the correct size
cd /private/groups/patenlab/mira/centrolign/analysis/HPRC_release2_QCv2/sample_lists/
ls  | while read line ; do wc -l $line;done
```

Copy scripts for plots to local
```sh
mkdir -p /private/groups/patenlab/mira/HPRC_release2_QCv2_all_pairs_heatmaps/
cp /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/*format5* /private/groups/patenlab/mira/HPRC_release2_QCv2_all_pairs_heatmaps/

cp /private/groups/patenlab/mira/centrolign/analysis/HPRC_release2_QCv2/sample_lists/*samples.txt /private/groups/patenlab/mira/HPRC_release2_QCv2_all_pairs_heatmaps/

cp /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/* /private/groups/patenlab/mira/HPRC_release2_QCv2_all_pairs_heatmaps/
```

Plot pairwise heatmaps
```sh
chromosomes=("chr13" "chr14" "chr17" "chr18")

for chr in "${chromosomes[@]}"
do
    python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap_v2.py \
        -t /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_all_pairs_nj_tree.format5.nwk \
        -s /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}.samples.txt \
        -p /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_pairwise_distance.csv \
        -m "Centrolign all pairs distances" \
        -n "${chr} NJ tree" \
        -d "All pairs Distances" \
        -o /Users/miramastoras/Desktop/github_repos/centrolign_analysis/analysis_notes/release2_QC_v2/figures/${chr}_r2_QC_v2_all_pairs --no_labels
done
```
## Dividing up the input guide trees for runtime

Because centrolign runtime gets unreasonable beyond ~ 180-200 samples, we need to divide up the input NJ trees and submit separate GSA runs. This works out for our downstream analysis okay though, because there are clades that have almost 0 alignment identity, we can't report variants between those clades, and will instead need to report variants only within clades that are alignable.

The goal right now is to loosely identify places we can break the tree, and safely keep samples apart that we are confident will not be alignable in the MSA.

Karen had an idea to also do a sanity check by aligning a smaller subset of samples between these clades and confirm there is not alignment between.

1. Adjustments to heatmap script, to pass in a list of samples and pain their pairwise combinations the same color (to visually check which samples got in the same grouping )

2. New script to take in the entire pairwise matrix, and a subset of samples, and print out for every sample IN the subset, the number of values it aligns to for samples OUTSIDE the subset, that have

### Starting with chr 5, 8, 12.

Chrom 5 is easy to begin with because it has a clear place to split the tree, at the deepest node.


Combine all fasta files for chr 5 first:
```sh
chromosomes=("chr5" "chr8" "chr12")

for chr in "${chromosomes[@]}"
do
  cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/fasta_lists/release2_QC_v2_${chr}.txt | while read line ; do cat $line ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.${chr}.fasta

  samtools faidx /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.${chr}.fasta

done

# Sanity check all seqs made it into the fasta by checking against sample list
for chr in "${chromosomes[@]}"
do
  list=`cat /private/groups/patenlab/mira/centrolign/analysis/HPRC_release2_QCv2/sample_lists/${chr}.samples.txt | wc -l`
  fasta=`cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.${chr}.fasta.fai | wc -l`

  echo $chr,$fasta,$list
  if [[ "$fasta" == "$list" ]]; then
    echo "true"
  fi
done

## Split into N subgroups
N=2
for chr in "${chromosomes[@]}"
do
  mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/${chr}/
  docker run -u `id -u`:`id -g` -v /private/groups:/private/groups/ \
      miramastoras/centromere_scripts:v0.1.4 \
      python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/split_fasta_by_tree.py \
      -t /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/${chr}_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
      -f /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.${chr}.fasta \
      -n ${N} \
      -o /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/${chr}/split_by_${N}/
done
```

Sanity check 1: color all samples in subgroup 0 with solid red
```sh
# get subgroup 0 list of samples for all three chroms
chromosomes=("chr5" "chr8" "chr12")

for chr in "${chromosomes[@]}"
do
  cut -f1 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/${chr}/split_by_2/subgroup_0_seqs.fasta.fai > /private/groups/patenlab/mira/color_subgroups_heatmap/${chr}.subgroup_0.samples.txt
done

```
Plotting on local computer:

```sh
chromosomes=("chr5" "chr8" "chr12")

for chr in "${chromosomes[@]}"
do
  python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap_v2.py \
    -t /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_all_pairs_nj_tree.format5.nwk \
    -s /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}.samples.txt \
    -p /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_pairwise_distance.csv \
    -m "Centrolign all pairs distances" \
    -n "${chr} NJ tree" \
    -d "All pairs Distances" \
    -o /Users/miramastoras/Desktop/color_subgroups_heatmap/${chr}_r2_QC_v2_all_pairs --no_labels \
    --highlight_samples /Users/miramastoras/Desktop/color_subgroups_heatmap/${chr}.subgroup_0.samples.txt
done
```

Sanity check 2: Bin pairwise distances between the two subset groups.
```sh
chromosomes=("chr5" "chr8" "chr12")
for chr in "${chromosomes[@]}"
do
    echo $chr
    python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/bin_pairwise_distances.py \
    -s /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}.samples.txt \
    -d /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_pairwise_distance.csv \
    -u /Users/miramastoras/Desktop/color_subgroups_heatmap/${chr}.subgroup_0.samples.txt \
    -o /Users/miramastoras/Desktop/github_repos/centrolign_analysis/analysis_notes/release2_QC_v2/split_nj_tree/${chr}_aligns_to_other_subgroups.csv
done
```
