## Simulating HOR sequences to test and evaluate centrolign


### 1. Identify cyclic HOR arrays without SVs

https://pmc.ncbi.nlm.nih.gov/articles/PMC9248890/

```
chr 2,3,4,6,7,10,11,12,14,15,16,17,20,21,22,X,Y
```

AS_HOR bed from Hailey for CHM13: https://public.gi.ucsc.edu/~hloucks/CenSat/CHM13/chm13v2.0.labels.as_hor.bed

Extract HOR from chm13 for each chromosome
```sh
chromosomes=("chr2" "chr3" "chr4" "chr6" "chr7" "chr10" "chr11" "chr12" "chr14" "chr15" "chr16" "chr17" "chr20" "chr21" "chr22" "chrX" "chrY")

mkdir -p per_chrom/work

for chr in "${chromosomes[@]}"
do
  # bedtools intersect alpha annotation with censat annotation to just get active hors
  echo "processing ${chr} "
  grep -w $chr chm13v2.0_censat_v2.1.bed | grep "hor" | grep -v "dhor" > per_chrom/work/chm13v2.0_censat_v2.1.${chr}.hor.bed

  bedtools intersect -wa -a chm13v2.0.labels.as_hor.bed -b per_chrom/work/chm13v2.0_censat_v2.1.${chr}.hor.bed > per_chrom/work/chm13v2.0.labels.as_hor.${chr}.active.bed

  # get start and end of active array from as_hor bed
  awk 'NR==1{start=$2; chrom=$1} END{print chrom"\t"start"\t"$3}' per_chrom/work/chm13v2.0.labels.as_hor.${chr}.active.bed > per_chrom/work/chm13v2.0.labels.as_hor.${chr}.active.mrg.bed

  # extract fasta sequence for active hor array
  bedtools getfasta -fi /private/groups/patenlab/mira/data/chm13v2.0.fa -bed per_chrom/work/chm13v2.0.labels.as_hor.${chr}.active.mrg.bed | sed '/^>/b; s/[a-z]/\U&/g' | sed '/^>/s/:.*//' > per_chrom/chm13v2.0.${chr}.active_hor.upper.fa

  # shift hor arrays over to start at 0
  awk 'NR==1 {shift=$2} {print $1"\t"($2-shift)"\t"($3-shift)"\t"$4}' per_chrom/work/chm13v2.0.labels.as_hor.${chr}.active.bed > per_chrom/chm13v2.0.labels.as_hor.${chr}.active.shifted.bed

done
```

### 2. Run MSA simulations


#### Step 1: make sim cases

```sh

mkdir -p /private/groups/patenlab/mira/centrolign/simulations/MSA_simulations/make_sim_cases_slurm_logs

git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/ pull

cd /private/groups/patenlab/mira/centrolign/simulations/MSA_simulations

chromosomes=("chr2" "chr3" "chr4" "chr6" "chr7" "chr10" "chr11" "chr12" "chr14" "chr15" "chr16" "chr17" "chr20" "chr21" "chr22" "chrX" "chrY")

for chr in "${chromosomes[@]}"
do
    sbatch /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/simulations/slurm_scripts/make_sim_cases_MSA.sh $chr
done
```

#### Step 2: Run centrolign MSA on simulated seqs and analyze results

```sh

mkdir -p /private/groups/patenlab/mira/centrolign/simulations/MSA_simulations/simulated_centrolign_slurm_logs

git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/ pull

cd /private/groups/patenlab/mira/centrolign/simulations/MSA_simulations

chromosomes=("chr2" "chr3" "chr4" "chr6" "chr7" "chr10" "chr11" "chr12" "chr14" "chr15" "chr16" "chr17" "chr20" "chr21" "chr22" "chrX" "chrY")

for chr in "${chromosomes[@]}"
do
    sbatch /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/simulations/slurm_scripts/simulated_centrolign_MSA.sh $chr
done
```
#### Plot results for MSA simulations

```sh
cd /private/groups/patenlab/mira/centrolign/simulations/MSA_simulations

chromosomes=("chr2" "chr3" "chr4" "chr6" "chr7" "chr10" "chr11" "chr12" "chr14" "chr15" "chr16" "chr17" "chr20" "chr21" "chr22" "chrX" "chrY")

mkdir -p /private/groups/patenlab/mira/centrolign/simulations/MSA_simulations/summary_tables

for chr in "${chromosomes[@]}"
do
    echo "Processing $chr"
    cat msa_${chr}_sim_cases_20250402/*/aln_summary_table.txt > summary_tables/msa_${chr}_sim_cases_20250402_aln_summary_tables.txt
done
```
Plot results
```sh
cd /private/groups/patenlab/mira/centrolign/simulations/MSA_simulations/summary_tables

ls | while read line ; do
    chr=`echo $line | cut -f2 -d"_"`
    Rscript /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/msa_simulations.R $line $chr /private/groups/patenlab/mira/centrolign/simulations/MSA_simulations/png_plots/${chr}_MSA_simulations
  done
```

## 3. Run pairwise simulations and comparisons to other tools

Location:
```
/private/groups/patenlab/mira/centrolign/simulations/pairwise_simulations
```

#### Step 1: make sim cases

```sh

mkdir -p /private/groups/patenlab/mira/centrolign/simulations/pairwise_simulations/make_sim_cases_slurm_logs

git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/ pull

cd /private/groups/patenlab/mira/centrolign/simulations/pairwise_simulations

chromosomes=("chr2" "chr3" "chr4" "chr6" "chr7" "chr10" "chr11" "chr12" "chr14" "chr15" "chr16" "chr17" "chr20" "chr21" "chr22" "chrX" "chrY")

for chr in "${chromosomes[@]}"
do
    sbatch /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/simulations/slurm_scripts/make_sim_cases_pairwise.sh $chr
done
```
#### Step 2: run centrolign and other pairwise aligners
```sh

mkdir -p /private/groups/patenlab/mira/centrolign/simulations/pairwise_simulations/simulated_centrolign_slurm_logs

chromosomes=("chr2" "chr3" "chr4" "chr6" "chr7" "chr10" "chr11" "chr12" "chr14" "chr15" "chr16" "chr17" "chr20" "chr21" "chr22" "chrX" "chrY")

git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/ pull

cd /private/groups/patenlab/mira/centrolign/simulations/pairwise_simulations

for chr in "${chromosomes[@]}"
do
    sbatch /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/simulations/slurm_scripts/simulated_centrolign_pairwise.sh $chr
done

# manually resubmitted # 60 for chr7
sbatch /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/simulations/slurm_scripts/simulated_centrolign_pairwise.sh chr7
```
/private/groups/patenlab/mira/centrolign/simulations/pairwise_simulations/pair_chr7_sim_cases_20250421/gen200/case_52

```
# manually resubmitted # 52 for chr14
sbatch /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/simulations/slurm_scripts/simulated_centrolign_pairwise.sh chr14

# 52 in array
pair_chr14_sim_cases_20250421/gen25/case_24/
```
#### 3. Plot results for pairwise simulations

```sh
cd /private/groups/patenlab/mira/centrolign/simulations/pairwise_simulations

chromosomes=("chr2" "chr3" "chr4" "chr6" "chr7" "chr10" "chr11" "chr12" "chr14" "chr15" "chr16" "chr17" "chr20" "chr21" "chr22" "chrX" "chrY")

mkdir -p /private/groups/patenlab/mira/centrolign/simulations/pairwise_simulations/summary_tables

for chr in "${chromosomes[@]}"
do
    echo "Processing $chr"
    cat pair_${chr}_sim_cases_20250421/*/*/aln_summary_table.txt > summary_tables/pair_${chr}_sim_cases_20250421_aln_summary_tables.txt
done
```
Plot results
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/simulations/pairwise_simulations/png_plots/

cd /private/groups/patenlab/mira/centrolign/simulations/pairwise_simulations/summary_tables

ls | while read line ; do
    chr=`echo $line | cut -f2 -d"_"`
    Rscript /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/pairwise_simulations.R $line $chr /private/groups/patenlab/mira/centrolign/simulations/pairwise_simulations/png_plots/${chr}_pairwise_simulations
  done
```

Investigating low chr11 performance
```
Rscript /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_simulations.R /Users/miramastoras/Desktop/pair_chr11_sim_cases_20250421_aln_summary_tables.txt chr11 /Users/miramastoras/Desktop/chr11_test

sort -k1,1 -k10,10 -k11,11 /Users/miramastoras/Desktop/pair_chr11_sim_cases_20250421_aln_summary_tables.txt > /Users/miramastoras/Desktop/pair_chr11_sim_cases_20250421_aln_summary_tables.sorted.txt
```
```
./plot_dotplot_alignment.py fasta1 fasta2 cigar[,cigar2,cigar3,...] svg_out_name
```
