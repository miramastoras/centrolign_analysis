### Using centrolign all pairs to create HOR trees as inputs

Generate NJ trees for all samples
```sh
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr7" "chr8" "chr9" "chr10" "chr11" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

for chr in "${chromosomes[@]}"
do
    docker run -u `id -u`:`id -g` -v /private/groups:/private/groups/ \
        miramastoras/centromere_scripts:v0.1.4 \
        python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/data_exploration/infer_tree.py \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/${chr}/pairwise_cigar/ \
        > /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_all_pairs/${chr}_r2_centrolign_all_pairs_nj_tree.nwk
  done 
```
