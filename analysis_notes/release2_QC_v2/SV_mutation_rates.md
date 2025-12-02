# Calculating mutation rates from centrolign SV calls

Mutation rate formula:

```sh
# SV calls in pairwise cigar / patristic distance } summed over all pairwise comparisons / # of pairwise comparisons
```
Patristic distances:

```sh
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

for chr in "${chromosomes[@]}"; do
  /private/home/mmastora/progs/centrolign/build/tree_pair_dist /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/${chr}_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk > /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/tree_pair_dist/${chr}_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk.pair_dists.tsv
done
```

Getting count of SVs per pairwise combination
```sh

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

for chr in "${chromosomes[@]}"; do
  for f in *.bed; do
    awk -F'\t' '
    NR==1 { s1=$1; s2=$4 }  # Capture sample1 and sample2 from first row
    {
        if ($8 == -1) c1++
        else if ($8 > 0 && $8 <= 0.10) c2++
        else if ($8 > 0.10) c3++
    }
    END {
        printf "%s\t%s\t%d\t%d\t%d\n", s1, s2, c1, c2, c3
    }' "$f"
  done > 
done
```
