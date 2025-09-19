## Comparing different formulas for converting a cigar string to distance

The purpose of this analysis is to find the best way to convert the cigar string to a pairwise distance for centrolign.

Formulas tested:

1. Jordan's original formula

```
1.0 - (2.0 * matches)
       ---------------
     (ref_len + query_len)
         |          |_ matches + mismatches + unaligned bases in query
         |_ matches + mismatches + unaligned bases in ref
```
2. Jordan's proposal involving the proportion of unaligned sequence:
```
identity:
(proportion aligned) * (matches / (matches + mismatches))

distance:
(proportion unaligned) * (mismatches / matches + mismatches)
```
3. Benedict's proposal:
```
mismatches/(matches+mismatches)
```
> We decided to scrap this one, because it breaks down in cases like this `210=2359579D3069392I168=` where the two sequences have a very small # of aligned bases

### 1. Implementing the formulas

Selecting two test cigar strings to operate on:

```
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/cigar_to_distance.py \
  -a /Users/miramastoras/Desktop/github_repos/centrolign_analysis/analysis_notes/cigar_to_dist/pairwise_cigar1/ \
  -o /Users/miramastoras/Desktop/github_repos/centrolign_analysis/analysis_notes/cigar_to_dist/formula2_ \
  -d 2

  python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/cigar_to_distance.py \
    -a /Users/miramastoras/Desktop/github_repos/centrolign_analysis/analysis_notes/cigar_to_dist/pairwise_cigar1/ \
    -o /Users/miramastoras/Desktop/github_repos/centrolign_analysis/analysis_notes/cigar_to_dist/formula4_ \
    -d 4
```

```
Cigar: 210=2359579D3069392I168=

Proportion aligned = (ref aligned + query aligned / ref total + query total)
1 - ((proportion aligned) * (matches / (mismatches + matches)))
1 - ( (  (378*2) / ((378*2)+2359579+3069392+ ) * (378 / (0 + 378))) = 0.999860766



Cigar: 915874=6X100=
Proportion aligned = (ref aligned + query aligned / ref total + query total)
Proportion aligned = 1
1 - ((proportion aligned) * (matches / (mismatches + matches)))
1 - (1 * (915974 / 915980 ) ) = 0.00000655036

Cigar: 100=2359D3060I1689=
Proportion aligned = (ref aligned + query aligned / ref total + query total)
Proportion aligned = ((100+1689)*2) / (((100+1689)*2) + 2359+3060) = 0.39768811826

1 - ((proportion aligned) * (matches / (mismatches + matches)))
1 - ((0.39768811826) * ((100+1689)/(100+1689))) = 0.60231188174

1.0 - (2.0 * matches) / (ref_len + query_len)
1 - (((100+1689)*2) / (((100+1689)*2) + 2359+3060)) = 0.60231188173

Cigar: 100=2359I16890=
Proportion aligned = (ref aligned + query aligned / ref total + query total)
Proportion aligned = ((16890+100) *2) / (((16890+100) *2) + 2359) = 0.93508351908

1 - ((proportion aligned) * (matches / (mismatches + matches)))
1 - ((0.93508351908) * (16990 / 16990)) = 0.06491648092
```
### 2. Plotting # mismatches vs alignment length

```
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/plot_cigar_dist.py -a /Users/miramastoras/Desktop/github_repos/centrolign_analysis/analysis_notes/cigar_to_dist/pairwise_cigar1/ -l test -o /Users/miramastoras/Desktop/github_repos/centrolign_analysis/analysis_notes/cigar_to_dist/test
```

Chr 3:
```
docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups \
  miramastoras/centromere_scripts:v0.1.4 python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/plot_cigar_dist.py -a /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr12/pairwise_cigar/ -l "HPRC release2 chr12" -o /private/groups/patenlab/mira/centrolign/distance_metric_testing/plot_cigar_dist/release2/chr12
```

Slurm script for all chromosomes:
```
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

for chr in "${chromosomes[@]}"
do
  echo $chr >> chr.txt
done
```


```sh
#!/bin/bash
# Slurm script to use centrolign as a pairwise aligner on all pairs of sequences
# from an input directory

#SBATCH --job-name=plot_cigar_dist
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=56gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=[1-24]%24
#SBATCH --output=logs/array_job_%A_task_%a.log
#SBATCH --time=1:00:00

CHR=$(awk "NR==$SLURM_ARRAY_TASK_ID" /private/groups/patenlab/mira/centrolign/distance_metric_testing/plot_cigar_dist/release2/chr.txt)

docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
  miramastoras/centromere_scripts:v0.1.4 \
  python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/plot_cigar_dist.py \
  -a /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/${CHR}/pairwise_cigar/ \
  -l "HPRC release2 ${CHR}" \
  -o /private/groups/patenlab/mira/centrolign/distance_metric_testing/plot_cigar_dist/release2/${CHR}
```
srun \
  --job-name "yak" \
  --cpus-per-task 4 \
  --partition short \
  --mem 50G \
  --time 1:00:00 \
  --pty bash \
  -i
