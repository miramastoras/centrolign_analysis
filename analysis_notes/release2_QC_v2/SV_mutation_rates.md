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
cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/centrolign_results.csv | grep "release 2 QC v2" | grep -v "chr20," | cut -f1 -d","   | while IFS=',' read -r subgroup ; do
  CHR=`echo $subgroup | cut -f1 -d"_"`
  echo $CHR, $subgroup

  cd /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/${CHR}/SV_beds/${subgroup}/

  for f in *.bed; do
    awk -F'\t' '
    NR==1 { s1=$1; s2=$4 }  # Capture sample1 and sample2 from first row
    {
        if ($8 == -1) c1++
      else if ($8 >= 0 && $8 <= 0.10) c2++
        else if ($8 > 0.10) c3++
    }
    END {
        printf "%s\t%s\t%d\t%d\t%d\n", s1, s2, c1, c2, c3
    }' "$f"
  done > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/${CHR}/${CHR}.${subgroup}.SV_summary.txt
done
```

Get list of these summary files
```sh
grep "release 2 QC v2" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/centrolign_results.csv | grep -v "chr20," | cut -f1 -d","   | while IFS=',' read -r subgroup ; do
  CHR=`echo $subgroup | cut -f1 -d"_"`

  echo /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/${CHR}/${CHR}.${subgroup}.SV_summary.txt
done > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/all_chrom_SV_summary.txt
```
Plotting in python notebook:

#### Plotting number of triangles inside tree heatmap

```sh
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

for chr in "${chromosomes[@]}"
do
  echo $chr
  grep $chr mutation_rates_scaled.csv | cut -f1,2,4 -d"," > mutation_rates_scaled.${chr}.triangles.csv
  grep $chr mutation_rates_scaled.csv | wc -l
done
```

```sh
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")


for chr in "${chromosomes[@]}"
do
    python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap_v2.py \
        -t /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_all_pairs_nj_tree.format5.nwk \
        -s /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}.samples.txt \
        -p /Users/miramastoras/Desktop/mutation_rates/mutation_rates_scaled.${chr}.triangles.csv \
        -m "Number of triangles" \
        -n "${chr} HOR NJ tree" \
        -d "All pairs Distances" \
        -o /Users/miramastoras/Desktop/mutation_rates/${chr}.triangles --no_labels
done
```
