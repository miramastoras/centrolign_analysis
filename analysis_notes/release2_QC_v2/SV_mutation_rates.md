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

## run it for cenhap trees that are complete
/private/home/mmastora/progs/centrolign/build/tree_pair_dist /private/groups/patenlab/mira/centrolign/annotations/hprc_cenhap_trees_12042025/chr11/full/HPRC_chr11_50650000_51023358_54476419_54808189_het449_m_final_dgp_upgma.nwk > /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/tree_pair_dist/cenhap_12042025/HPRC_chr11_50650000_51023358_54476419_54808189_het449_m_final_dgp_upgma.pair_dists.tsv

/private/home/mmastora/progs/centrolign/build/tree_pair_dist /private/groups/patenlab/mira/centrolign/annotations/hprc_cenhap_trees_12042025/chr12_prelim/HPRC_chr12_34544731_34593492_37202490_37285321_het60_m_hprc_dgp_rnj_upgma.nwk > /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/tree_pair_dist/cenhap_12042025/HPRC_chr12_34544731_34593492_37202490_37285321_het60_m_hprc_dgp_rnj_upgma.pair_dists.tsv

/private/home/mmastora/progs/centrolign/build/tree_pair_dist /private/groups/patenlab/mira/centrolign/annotations/hprc_cenhap_trees_12042025/chr17_prelim/HPRC_chr17_23278614_23433372_27571319_27700000_het114_m_hprc_dgp_upgma.nwk > /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/tree_pair_dist/cenhap_12042025/HPRC_chr17_23278614_23433372_27571319_27700000_het114_m_hprc_dgp_upgma.pair_dists.tsv

/private/home/mmastora/progs/centrolign/build/tree_pair_dist /private/groups/patenlab/mira/centrolign/annotations/hprc_cenhap_trees_12042025/chr6_prelim/standard/HPRC_chr6_58200000_58286706_61058390_61123742_het70_m_hprc_dgp_rnj_upgma.nwk > /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/tree_pair_dist/cenhap_12042025/HPRC_chr6_58200000_58286706_61058390_61123742_het70_m_hprc_dgp_rnj_upgma.pair_dists.tsv
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

Plotting normalized number of triangles

```sh
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

for chr in "${chromosomes[@]}"
do
  echo $chr
  grep $chr, /private/groups/patenlab/mira/result_split_samples.csv | grep "triangle" | cut -f1,2,3 -d"," > results_split_samples/mutation_rates_scaled.${chr}.triangles.csv
  grep $chr, /private/groups/patenlab/mira/result_split_samples.csv | wc -l
done
```

```sh
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")


for chr in "${chromosomes[@]}"
do
    python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap_v2.py \
        -t /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_all_pairs_nj_tree.format5.nwk \
        -s /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}.samples.txt \
        -p /Users/miramastoras/Desktop/results_split_samples/mutation_rates_scaled.${chr}.triangles.csv \
        -m "Number of triangles" \
        -n "${chr} HOR NJ tree" \
        -d "All pairs Distances" \
        -o /Users/miramastoras/Desktop/results_split_samples/${chr}.triangles --no_labels
done
```

#### Testing using a distance thats just based off the SNPs in the array - this will only work for very short distances

Run cigar to dist conversion for all chromosomes

```sh
#!/bin/bash
#SBATCH --job-name=cigar_to_dist
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --array=[0-23]%24
#SBATCH --time=12:00:00

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

chr=${chromosomes[$SLURM_ARRAY_TASK_ID]}

python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/cigar_to_distance.py -d 3 \
    -a /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/${chr}/pairwise_cigar/ \
    -o /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/SNP_based_distances/${chr}_r2_QC_v2_centrolign_
```

Submitting script:
```sh
cd /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/SNP_based_distances
mkdir -p logs
sbatch slurm_cigar_to_dist.sh
```
