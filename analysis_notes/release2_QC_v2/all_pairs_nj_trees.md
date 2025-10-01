## Create neighbor-joining trees based on all pairs alignments, and decide where to split trees for MSA runs

### Chr 1, Chr 12, Chr 10, Chr 11, Chr 8

Get distance matrices for all pairs
```sh

## Submit on slurm - ~ 1hr runtime
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
