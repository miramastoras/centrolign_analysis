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

Slurm script to run on list of samples, downloads fasta file per sample and extracts HOR sequence per sample, placing it in dir per chromosome
```sh

git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas

mkdir -p logs
sbatch \
  /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/extract_fasta_r2_QCv2.sh
```
Accounting check: Make sure we aren't missing any fastas
```sh
# check to make sure none of the fastas are empty
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas

find . -type f -empty

# check that number of fastas per chrom matches julian's original files
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

for chr in "${chromosomes[@]}"
do
  fastas=`ls /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/${chr}/ | wc -l`
  QC=`grep -v "sample_id" /private/groups/migalab/juklucas/censat_regions/active_arrays/asat_arrays_${chr}.tsv | wc -l`

  echo $chr,$fastas,$QC
  if [[ "$fastas" == "$QC" ]]; then
    echo "true"
  fi
done

## checks passed.
```
### Run centrolign all pairs for Chr 1, Chr 12, Chr 10, Chr 11, Chr 8


Get list of fasta files per chromosome
```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/fasta_lists/

chromosomes=("chr1" "chr12" "chr10" "chr11" "chr8")
for chr in "${chromosomes[@]}"
do
    ls /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/${chr}/ | while read line ;
      do realpath $chr/$line >> /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/fasta_lists/release2_QC_v2_${chr}.txt
    done
done

# Sanity check again that we have all fastas
for chr in "${chromosomes[@]}"
do
  fastas=`cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/fasta_lists/release2_QC_v2_${chr}.txt | wc -l`
  QC=`grep -v "sample_id" /private/groups/migalab/juklucas/censat_regions/active_arrays/asat_arrays_${chr}.tsv | wc -l`

  echo $chr,$fastas,$QC
  if [[ "$fastas" == "$QC" ]]; then
    echo "true"
  fi
done
```
Generate all vs all combinations
```sh
chromosomes=("chr1" "chr12" "chr10" "chr11" "chr8")

# get all vs all pairwise combinations
for chr in "${chromosomes[@]}"
do
  mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/${chr}

  while read -r s1; do
      while read -r s2; do
          [[ "$s1" < "$s2" ]] && echo -e "$s1\t$s2\t$chr"
      done < /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/fasta_lists/release2_QC_v2_${chr}.txt
  done < /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/fasta_lists/release2_QC_v2_${chr}.txt > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/${chr}/release2_QC_v2_all_pairs_combinations_${chr}.txt
done
```
Count per chrom
```sh
for chr in "${chromosomes[@]}"
do
  wc -l /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/${chr}/release2_QC_v2_all_pairs_combinations_${chr}.txt
done

# 63903 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr1/release2_QC_v2_all_pairs_combinations_chr1.txt
# 69378 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr12/release2_QC_v2_all_pairs_combinations_chr12.txt
# 57630 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr10/release2_QC_v2_all_pairs_combinations_chr10.txt
# 84666 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr11/release2_QC_v2_all_pairs_combinations_chr11.txt
# 61075 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr8/release2_QC_v2_all_pairs_combinations_chr8.txt
```

Run centrolign all pairs - Chr 1
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull


mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr1/logs

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr1/

sbatch \
    --job-name=chr1_r2_QCv2 \
    --array=[1-30000]%128 \
    --export=CHR=chr1 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh

sbatch \
    --job-name=chr1_r2_QCv2 \
    --array=[30001-63903]%128 \
    --export=CHR=chr1 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh
```

Run centrolign all pairs - Chr 12
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr12/logs

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr12/

sbatch \
    --job-name=chr12_r2_QCv2 \
    --array=[1-30000]%128 \
    --export=CHR=chr12 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh

sbatch \
    --job-name=chr12_r2_QCv2 \
    --array=[30001-69378]%128 \
    --export=CHR=chr12 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh
```
Run centrolign all pairs - Chr 10
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr10/logs

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr10/

sbatch \
    --job-name=chr10_r2_QCv2 \
    --array=[1-30000]%128 \
    --export=CHR=chr10 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh

# not submitted yet 
sbatch \
    --job-name=chr10_r2_QCv2 \
    --array=[30001-57630]%128 \
    --export=CHR=chr10 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh
```

**Sanity check: search for empty cigar strings, count number in each dir**
