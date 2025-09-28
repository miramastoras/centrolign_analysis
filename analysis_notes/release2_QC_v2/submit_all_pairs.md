## Submitting all vs all centrolign direct pairwise alignments for release2 QC v2

"release2 QCv2" indicates samples passing Julian's QC pipeline, which integrates Flagger and
nucFlag to filter arrays.


Downloaded full assembly list for release 2 from here https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/data_tables/assemblies_pre_release_v0.6.1.index.csv

```sh
/private/groups/patenlab/mira/centrolign/annotations/assemblies_pre_release_v0.6.1.index.csv
```
Manually added CHM13 and HG002 to the file

#### 1. Extract HOR fasta files for all chromosomes

Concatenate all of the QC csv files together. Have to run chr3 and chr4 separately because they are formatted differently

```sh
chromosomes=("chr1" "chr2" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

grep sample_id /private/groups/migalab/juklucas/censat_regions/active_arrays/asat_arrays_chr1.csv > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/asat_arrays_wo_chr3_4.csv

for chr in "${chromosomes[@]}"
do
  wc -l /private/groups/migalab/juklucas/censat_regions/active_arrays/asat_arrays_${chr}.csv
  cat /private/groups/migalab/juklucas/censat_regions/active_arrays/asat_arrays_${chr}.csv | grep -v "sample_id" >> /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/asat_arrays_wo_chr3_4.csv
done  

#!/bin/sh
chromosomes=("chr3" "chr4")

grep sample_id /private/groups/migalab/juklucas/censat_regions/active_arrays/asat_arrays_chr3.csv > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/asat_arrays_chr3_4.csv

for chr in "${chromosomes[@]}"
do
  wc -l /private/groups/migalab/juklucas/censat_regions/active_arrays/asat_arrays_${chr}.csv
  cat /private/groups/migalab/juklucas/censat_regions/active_arrays/asat_arrays_${chr}.csv | grep -v "sample_id" >> /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/asat_arrays_chr3_4.csv
done  
```

Python script to create per sample bed files containing all of the arrays
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/parse_QC_csv.py \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/asat_arrays_wo_chr3_4.csv \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/

# Chr 3 and Chr 4
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/parse_QC_csv.py \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/asat_arrays_chr3_4.csv \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/
```

Get list of all samples for submitting sbatch script
```sh
cut -f1-3 -d"," asat_arrays_chr3_4.csv | grep -v sample_id > all_samples_with_asats.txt
cut -f1-3 -d"," asat_arrays_wo_chr3_4.csv | grep -v sample_id >> all_samples_with_asats.txt

sort all_samples_with_asats.txt | uniq > tmp ; mv tmp all_samples_with_asats.txt
# 466 samples
```

Slurm script to run on list of samples, downloads fasta file per sample and extracts HOR sequence per sample placing it in dir per chromosome
```sh
sbatch \
  extract_fasta_r2_QCv2.sh

```

Get list of fasta file per chromosome, and all vs all combinations









```sh
# Python script to create new csv file with s3 links, and output a bed file
# with asat locations

git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

chromosomes=("chr1" "chr12" "chr10" "chr11" "chr8")
for chr in "${chromosomes[@]}"
do
    mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/active_array_fastas/${chr}/

    mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/active_array_beds/${chr}/

    python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/parse_QC_csv.py \
      /private/groups/migalab/juklucas/censat_regions/active_arrays/asat_arrays_${chr}.csv \
      /private/groups/patenlab/mira/centrolign/annotations/assemblies_pre_release_v0.6.1.index.csv \
      /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/active_array_csvs/asat_arrays_${chr}.csv \
      /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/active_array_beds/${chr}/
done

# Cat all of the active array files together
chromosomes=("chr1" "chr12" "chr10" "chr11" "chr8")
for chr in "${chromosomes[@]}"
do
  cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/active_array_csvs/asat_arrays_${chr}.csv >> /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/active_array_csvs/asat_arrays_all_chroms.csv
done
```
Run slurm script to extract fasta file


Need to use python script to get list of all arrays passing QC per sample. Just extract fastas for all chroms, then just run the all pairs for the first 5.


### Starting with Chr 1, Chr 12, Chr 10, Chr 11, Chr 8
