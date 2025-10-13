## Submitting all vs all centrolign direct pairwise alignments for release2 QC v2

"release2 QCv2" indicates samples passing Julian's QC pipeline, which integrates Flagger and
nucFlag to filter arrays.


Downloaded full assembly list for release 2 from here https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/data_tables/assemblies_pre_release_v0.6.1.index.csv

```sh
/private/groups/patenlab/mira/centrolign/annotations/assemblies_pre_release_v0.6.1.index.csv
```
Manually added CHM13 and HG002 to the file

#### 1. Extract HOR fasta files for all chromosomes

Concatenate all of the QC csv files together. Have to run chr3 and chr4 separately because they came in later.

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

grep sample_id /private/groups/migalab/juklucas/censat_regions/active_arrays/asat_arrays_chr3.tsv > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/asat_arrays_chr3_4.tsv

for chr in "${chromosomes[@]}"
do
  wc -l /private/groups/migalab/juklucas/censat_regions/active_arrays/asat_arrays_${chr}.tsv
  cat /private/groups/migalab/juklucas/censat_regions/active_arrays/asat_arrays_${chr}.tsv | grep -v "sample_id" >> /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/asat_arrays_chr3_4.tsv
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
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/asat_arrays_chr3_4.tsv \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/
```

Get list of all samples for submitting sbatch script
```sh
cut -f1-3 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/asat_arrays_chr3_4.tsv | grep -v sample_id | sed 's/\t/,/g' > all_samples_with_asats.txt
cut -f1-3 -d"," asat_arrays_wo_chr3_4.csv | grep -v sample_id >> all_samples_with_asats.txt

sort all_samples_with_asats.txt | uniq > tmp ; mv tmp all_samples_with_asats.txt
# 465 samples
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

sbatch \
    --job-name=chr10_r2_QCv2 \
    --array=[30001-57630]%128 \
    --export=CHR=chr10 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh
```

Run centrolign all pairs - Chr 11
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr11/logs

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr11/

sbatch \
    --job-name=chr11_r2_QCv2 \
    --array=[1-40000]%128 \
    --export=CHR=chr11 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh

sbatch \
    --job-name=chr11_r2_QCv2 \
    --array=[40001-84666]%128 \
    --export=CHR=chr11 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh
```
Run centrolign all pairs - Chr 8
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr8/logs

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr8/

sbatch \
    --job-name=chr8_r2_QCv2 \
    --array=[1-30000]%128 \
    --export=CHR=chr8 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh

sbatch \
    --job-name=chr8_r2_QCv2 \
    --array=[30001-61075]%128 \
    --export=CHR=chr8 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh
```

**Sanity check: search for empty cigar strings, count number in each dir**

```sh
# Search for empty cigar strings
chromosomes=("chr1" "chr12" "chr10" "chr11" "chr8")

for chr in "${chromosomes[@]}"
do
  echo $chr
  find /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/${chr}/pairwise_cigar/ -type f -empty
done

# Count strings in each dir
chromosomes=("chr1" "chr12" "chr10" "chr11" "chr8")

for chr in "${chromosomes[@]}"
do
  echo $chr
  ls /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/${chr}/pairwise_cigar/ | wc -l
done
```

### Run centrolign all pairs for Chr 2, Chr 5, Chr 6, Chr 7, Chr 9


Get list of fasta files per chromosome
```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/fasta_lists/

chromosomes=("chr2" "chr5" "chr6" "chr7" "chr9")
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
chromosomes=("chr2" "chr5" "chr6" "chr7" "chr9")

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

# 75078 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr2/release2_QC_v2_all_pairs_combinations_chr2.txt
# 18721 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr6/release2_QC_v2_all_pairs_combinations_chr6.txt
# 43956 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr7/release2_QC_v2_all_pairs_combinations_chr7.txt
# 56953 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr5/release2_QC_v2_all_pairs_combinations_chr5.txt
# 83028 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr9/release2_QC_v2_all_pairs_combinations_chr9.txt
```
Run centrolign all pairs - Chr 2
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr2/logs

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr2/

sbatch \
    --job-name=chr2_r2_QCv2 \
    --array=[1-40000]%128 \
    --export=CHR=chr2 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh

sbatch \
    --job-name=chr2_r2_QCv2 \
    --array=[40001-75078]%128 \
    --export=CHR=chr2 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh
```
Run centrolign all pairs - Chr 6
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr6/logs

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr6/

sbatch \
    --job-name=chr6_r2_QCv2 \
    --array=[1-18721]%128 \
    --export=CHR=chr6 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh
```

Run centrolign all pairs - Chr 7
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr7/logs

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr7/

sbatch \
    --job-name=chr7_r2_QCv2 \
    --array=[1-43956]%128 \
    --export=CHR=chr7 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh
```
Run centrolign all pairs - Chr 5
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr5/logs

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr5/

sbatch \
    --job-name=chr5_r2_QCv2 \
    --array=[1-43956]%128 \
    --export=CHR=chr5 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh

sbatch \
--job-name=chr5_r2_QCv2 \
    --array=[43957-56953]%128 \
    --export=CHR=chr5 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh
```
Run centrolign all pairs - Chr 9
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr9/logs

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr9/

sbatch \
    --job-name=chr9_r2_QCv2 \
    --array=[1-40000]%128 \
    --export=CHR=chr9 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh

#!/bin/sh
sbatch \
    --job-name=chr9_r2_QCv2 \
    --array=[40001-83028]%128 \
    --export=CHR=chr9 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh
```

**Sanity check: search for empty cigar strings, count number in each dir**

```sh
# Search for empty cigar strings
chromosomes=("chr2" "chr6" "chr7" "chr5" "chr9")

for chr in "${chromosomes[@]}"
do
  echo $chr
  find /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/${chr}/pairwise_cigar/ -type f -empty
done

# Count strings in each dir
chromosomes=("chr2" "chr6" "chr7" "chr5" "chr9")

for chr in "${chromosomes[@]}"
do
  echo $chr
  ls /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/${chr}/pairwise_cigar/ | wc -l
done
```

### Run centrolign all pairs for Chr 13, Chr 14, Chr 15, Chr 16, Chr 17, CHR 18


Get list of fasta files per chromosome
```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/fasta_lists/

chromosomes=("chr13" "chr14" "chr15" "chr16" "chr17" "chr18")
for chr in "${chromosomes[@]}"
do
    ls /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/${chr}/ | while read line ;
      do realpath $chr/$line >> /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/fasta_lists/release2_QC_v2_${chr}.txt
    done
done

# Sanity check again that we have all fastas - check against QC list length

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
chromosomes=("chr13" "chr14" "chr15" "chr16" "chr17" "chr18")

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

# 27730 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr13/release2_QC_v2_all_pairs_combinations_chr13.txt
# 72010 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr14/release2_QC_v2_all_pairs_combinations_chr14.txt
# 52003 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr15/release2_QC_v2_all_pairs_combinations_chr15.txt
# 76245 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr16/release2_QC_v2_all_pairs_combinations_chr16.txt
# 9591 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr17/release2_QC_v2_all_pairs_combinations_chr17.txt
# 6786 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr18/release2_QC_v2_all_pairs_combinations_chr18.txt
```
Run centrolign all pairs - Chr 13
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr13/logs

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr13/

sbatch \
    --job-name=chr13_r2_QCv2 \
    --array=[1-27730]%128 \
    --export=CHR=chr13 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh
```
Run centrolign all pairs - Chr 14
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr14/logs

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr14/

sbatch \
    --job-name=chr14_r2_QCv2 \
    --array=[1-30000]%128 \
    --export=CHR=chr14 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh

sbatch \
    --job-name=chr14_r2_QCv2 \
    --array=[30001-72010]%128 \
    --export=CHR=chr14 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh
```

Run centrolign all pairs - Chr 17
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr17/logs

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr17/

sbatch \
    --job-name=chr17_r2_QCv2 \
    --array=[1-9591]%128 \
    --export=CHR=chr17 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh
```

Run centrolign all pairs - Chr 18
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr18/logs

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr18/

sbatch \
    --job-name=chr18_r2_QCv2 \
    --array=[1-6786]%128 \
    --export=CHR=chr18 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh
```

Run centrolign all pairs - Chr 15
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr15/logs

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr15/

sbatch \
    --job-name=chr15_r2_QCv2 \
    --array=[1-52003]%128 \
    --export=CHR=chr15 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh
```
Run centrolign all pairs - Chr 16
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr16/logs

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr16/

sbatch \
    --job-name=chr16_r2_QCv2 \
    --array=[1-30000]%128 \
    --export=CHR=chr16 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh

sbatch \
    --job-name=chr16_r2_QCv2 \
    --array=[30001-76245]%128 \
    --export=CHR=chr16 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh
```

**Sanity check: search for empty cigar strings, count number in each dir**

```sh
# Search for empty cigar strings
chromosomes=("chr13" "chr14" "chr17" "chr18" "chr16" "chr15")

for chr in "${chromosomes[@]}"
do
  echo $chr
  find /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/${chr}/pairwise_cigar/ -type f -empty
done

# Count strings in each dir
chromosomes=("chr13" "chr14" "chr17" "chr18" "chr16" "chr15")

for chr in "${chromosomes[@]}"
do
  echo $chr
  ls /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/${chr}/pairwise_cigar/ | wc -l
done
```

### Run centrolign all pairs for Chr 19, Chr 20, Chr 21, Chr 22, Chr X, Chr Y


Get list of fasta files per chromosome
```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/fasta_lists/

chromosomes=("chr19" "chr20" "chr21" "chr22" "chrX" "chrY")
for chr in "${chromosomes[@]}"
do
    ls /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/${chr}/ | while read line ;
      do realpath $chr/$line >> /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/fasta_lists/release2_QC_v2_${chr}.txt
    done
done

# Sanity check again that we have all fastas - check against QC list length

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
chromosomes=("chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

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

# 48828 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr19/release2_QC_v2_all_pairs_combinations_chr19.txt
# 53956 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr20/release2_QC_v2_all_pairs_combinations_chr20.txt
# 38226 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr21/release2_QC_v2_all_pairs_combinations_chr21.txt
# 67161 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr22/release2_QC_v2_all_pairs_combinations_chr22.txt
# 47278 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chrX/release2_QC_v2_all_pairs_combinations_chrX.txt
# 990 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chrY/release2_QC_v2_all_pairs_combinations_chrY.txt
```
Run centrolign all pairs - Chr 19
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr19/logs

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr19/

sbatch \
    --job-name=chr19_r2_QCv2 \
    --array=[1-48828]%128 \
    --export=CHR=chr19 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh
```

Run centrolign all pairs - Chr 20
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr20/logs

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr20/

sbatch \
    --job-name=chr20_r2_QCv2 \
    --array=[1-30000]%128 \
    --export=CHR=chr20 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh

sbatch \
    --job-name=chr20_r2_QCv2 \
    --array=[30001-53956]%128 \
    --export=CHR=chr20 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh
```

Run centrolign all pairs - Chr 21
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr21/logs

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr21/

sbatch \
    --job-name=chr21_r2_QCv2 \
    --array=[1-38226]%128 \
    --export=CHR=chr21 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh
```

Run centrolign all pairs - Chr 22
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr22/logs

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr22/

sbatch \
    --job-name=chr22_r2_QCv2 \
    --array=[1-30000]%128 \
    --export=CHR=chr22 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh

sbatch \
    --job-name=chr22_r2_QCv2 \
    --array=[30001-67161]%128 \
    --export=CHR=chr22 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh
```
Run centrolign all pairs - Chr X
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chrX/logs

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chrX/

sbatch \
    --job-name=chrX_r2_QCv2 \
    --array=[1-47278]%128 \
    --export=CHR=chrX \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh
```
Run centrolign all pairs - Chr Y
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chrY/logs

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chrY/

sbatch \
    --job-name=chrY_r2_QCv2 \
    --array=[1-990]%128 \
    --export=CHR=chrY \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs.sh
```
**Sanity check: search for empty cigar strings, count number in each dir**

```sh
# Search for empty cigar strings
chromosomes=("chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

for chr in "${chromosomes[@]}"
do
  echo $chr
  find /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/${chr}/pairwise_cigar/ -type f -empty
done

# Count strings in each dir
chromosomes=("chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

for chr in "${chromosomes[@]}"
do
  echo $chr
  ls /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/${chr}/pairwise_cigar/ | wc -l
done
```
