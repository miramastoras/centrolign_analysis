# Preparing inputs to centrolign

Aligned to CHM13, then extracted HOR arrays that are contiguous
https://docs.google.com/spreadsheets/d/1is_jiWsDoqj_1QIcunGJLoojmToX3z9TFvTSEWCY2jA/edit?gid=1520027858#gid=1520027858

https://github.com/miramastoras/centrolign_analysis/tree/main/batch_submissions/Asm2AsmAlignerPaf/release2
https://github.com/miramastoras/centrolign_analysis/tree/main/batch_submissions/extract_hors_HPRC/release2

For each chr, get list of samples with contiguous arrays
```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/contiguous_HORs/

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

for chr in "${chromosomes[@]}"
do
    echo "Processing $chr"
    cut -f 1 -d"," extract_hors_HPRC_release2.csv | grep -v "sample_id" | \
    while read line ; do
        ls $line/analysis/extract_hors_HPRC_outputs/*${chr}_hor_array.fasta | cut -f4 -d"/" | cut -f3 -d"_"
      done | grep -v "cannot access" > /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/contiguous_HORs/HPRC_release2_contiguous_HORs_${chr}.txt
done

# Count number of contiguous arrays per chromosome
for chr in "${chromosomes[@]}"
do  
  count=`wc -l contiguous_HORs/HPRC_release2_contiguous_HORs_${chr}.txt | cut -f1 -d" "`
  echo $chr $count
done
```
```
chr1 398
chr2 433
chr3 160
chr4 196
chr5 375
chr6 191
chr7 337
chr8 358
chr9 407
chr10 356
chr11 422
chr12 385
chr13 303
chr14 348
chr15 331
chr16 402
chr17 139
chr18 130
chr19 392
chr20 343
chr21 275
chr22 333
chrX 315
chrY 49
```
Parse logs for reason for filtering
```
cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/slurm_logs

for log in *.log ;
  do grep "filtering" $log | grep -v "\["
done > filter_log.txt

# 491 filtering for gaps
# 2013 	filtering for discontiguous array
# 7 filtering for discordant mapped chromosome chr18
# 199 filtering HOR group for incomplete assignments
# 809 filtering for discontiguous array over chromosome _ in a resolved HOR group
# 10 filtering for being at end of contig of length _ in unique type
# 14 filtering for being at end of contig of length _ in a resolved HOR group

grep -v "filtering for discontiguous array over chromosome" filter_log.txt | grep -v "filtering for being at end of contig of length" | grep -v "filtering HOR group for incomplete assignments" | grep -v "filtering for discordant mapped chromosome" | grep -v "filtering for discontiguous array" | grep -v "filtering for gaps"

for log in *.log ;
  do grep "filtering for discordant HOR type" $log | grep -v "\["
done > filter_log2.txt
```
Get list of fasta paths per chromosome
```sh

cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

# get file containing all the filepaths to the fastas for each chr
for chr in "${chromosomes[@]}"
do
    echo "Processing $chr"
    cut -f 1 -d"," extract_hors_HPRC_release2.csv | grep -v "sample_id" | \
    while read line ; do
        ls $line/analysis/extract_hors_HPRC_outputs/*${chr}_hor_array.fasta
      done | grep -v "cannot access" > /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/HPRC_release2_contiguous_HORs_${chr}.fasta_list.txt
done

# combine all HOR sequence into one fasta per chromosome
cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2

for chr in "${chromosomes[@]}"
do
    cat /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/HPRC_release2_contiguous_HORs_${chr}.fasta_list.txt | while read line ; do cat $line ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/centrolign_fastas/HPRC_release2_HORs_${chr}.fasta
    samtools faidx /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/centrolign_fastas/HPRC_release2_HORs_${chr}.fasta
done
```
## Submit pairwise alignments

### Starting with Chr1, Chr6, Chr8, Chr10, Chr12

```sh
chromosomes=("chr1" "chr6" "chr8" "chr10" "chr12")

# make directory for analysis
for chr in "${chromosomes[@]}"
do
    mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/${chr}
    mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/${chr}/logs/
done

for chr in "${chromosomes[@]}"
do
  while read -r s1; do
      while read -r s2; do
          [[ "$s1" < "$s2" ]] && echo -e "$s1\t$s2"
      done < /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/contiguous_HORs/HPRC_release2_contiguous_HORs_${chr}.txt
  done < /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/contiguous_HORs/HPRC_release2_contiguous_HORs_${chr}.txt > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/${chr}/HPRC_release2_contiguous_HOR_all_pairs_combinations_${chr}.txt
done
```
Get number of combinations for each chrom
```sh
for chr in "${chromosomes[@]}"
do
  wc -l /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/${chr}/HPRC_release2_contiguous_HOR_all_pairs_combinations_${chr}.txt
done
```
```
79003 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr1/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr1.txt
18145 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr6/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr6.txt
63903 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr8/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr8.txt
63190 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr10/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr10.txt
73920 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr12/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr12.txt
```
Slurm parameters:
```
#SBATCH --job-name=chr6_pairwise-centrolign
#SBATCH --array=[1-18145]%128
CHR=chr6

#SBATCH --job-name=chr1_pairwise-centrolign
#SBATCH --array=[1-30000]%128
#SBATCH --array=[30001-78913]%128
#SBATCH --array=[78913-79003]%128
CHR=chr1

#SBATCH --job-name=chr12_pairwise-centrolign
#SBATCH --array=[1-30000]%128
#SBATCH --array=[30001-73829]%128
#SBATCH --array=[73829-73920]%128
CHR=chr12

#SBATCH --job-name=chr10_pairwise-centrolign
#SBATCH --array=[1-30000]%128
#SBATCH --array=[30001-63118]%128
#SBATCH --array=[63118-63190]%128
CHR=chr10

#SBATCH --job-name=chr8_pairwise-centrolign
#SBATCH --array=[1-30000]%128
#SBATCH --array=[30001-63823]%128
#SBATCH --array=[63823-63903]%128
CHR=chr8
```

### Starting with Chr2, Chr3, Chr4, Chr5, Chr7

```sh
chromosomes=("chr2" "chr3" "chr4" "chr5" "chr7")

# make directory for analysis
for chr in "${chromosomes[@]}"
do
    mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/${chr}
    mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/${chr}/logs/
done

for chr in "${chromosomes[@]}"
do
  while read -r s1; do
      while read -r s2; do
          [[ "$s1" < "$s2" ]] && echo -e "$s1\t$s2"
      done < /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/contiguous_HORs/HPRC_release2_contiguous_HORs_${chr}.txt
  done < /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/contiguous_HORs/HPRC_release2_contiguous_HORs_${chr}.txt > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/${chr}/HPRC_release2_contiguous_HOR_all_pairs_combinations_${chr}.txt
done
```
Get number of combinations for each chrom
```sh
for chr in "${chromosomes[@]}"
do
  wc -l /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/${chr}/HPRC_release2_contiguous_HOR_all_pairs_combinations_${chr}.txt
done
```
```
93528 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr2/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr2.txt
12720 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr3/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr3.txt
19110 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr4/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr4.txt
70125 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr5/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr5.txt
56616 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr7/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr7.txt
```
Slurm parameters:
```
#SBATCH --job-name=chr3_pairwise-centrolign
#SBATCH --array=[1-12720]%128
CHR=chr3

#SBATCH --job-name=chr4_pairwise-centrolign
#SBATCH --array=[1-19110]%128
CHR=chr4

#SBATCH --job-name=chr7_pairwise-centrolign
#SBATCH --array=[1-56616]%128
CHR=chr7

#SBATCH --job-name=chr5_pairwise-centrolign
#SBATCH --array=[1-30000]%128
#SBATCH --array=[30001-70125]%128
CHR=chr5

#SBATCH --job-name=chr2_pairwise-centrolign
#SBATCH --array=[1-30000]%128
#SBATCH --array=[30001-60000]%128
#SBATCH --array=[60001-93528]%128
CHR=chr2
```

### chr9" "chr11" "chr13" "chr14" "chr15

```sh
chromosomes=("chr9" "chr11" "chr13" "chr14" "chr15")

# make directory for analysis
for chr in "${chromosomes[@]}"
do
    mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/${chr}
    mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/${chr}/logs/
done

for chr in "${chromosomes[@]}"
do
  while read -r s1; do
      while read -r s2; do
          [[ "$s1" < "$s2" ]] && echo -e "$s1\t$s2"
      done < /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/contiguous_HORs/HPRC_release2_contiguous_HORs_${chr}.txt
  done < /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/contiguous_HORs/HPRC_release2_contiguous_HORs_${chr}.txt > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/${chr}/HPRC_release2_contiguous_HOR_all_pairs_combinations_${chr}.txt
done
```
Get number of combinations for each chrom
```sh
for chr in "${chromosomes[@]}"
do
  wc -l /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/${chr}/HPRC_release2_contiguous_HOR_all_pairs_combinations_${chr}.txt
done
```
```
82621 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr9/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr9.txt
88831 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr11/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr11.txt
45753 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr13/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr13.txt
60378 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr14/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr14.txt
54615 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr15/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr15.txt
```
Slurm parameters:
```
#SBATCH --job-name=chr9_pairwise-centrolign
#SBATCH --array=[1-40000]%128
#SBATCH --array=[40001-82621]%128
CHR=chr9

#SBATCH --job-name=chr11_pairwise-centrolign
#SBATCH --array=[1-40000]%128
#SBATCH --array=[40001-88831]%128
CHR=chr11

#SBATCH --job-name=chr13_pairwise-centrolign
#SBATCH --array=[1-45753]%128
CHR=chr13

#SBATCH --job-name=chr14_pairwise-centrolign
#SBATCH --array=[1-30000]%128
#SBATCH --array=[30001-60378]%128
CHR=chr14

#SBATCH --job-name=chr15_pairwise-centrolign
#SBATCH --array=[1-30000]%128
#SBATCH --array=[30001-54615]%128
CHR=chr15
```

### chr16" "chr17" "chr13" "chr14" "chr15

```sh
chromosomes=("chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

# make directory for analysis
for chr in "${chromosomes[@]}"
do
    mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/${chr}
    mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/${chr}/logs/
done

for chr in "${chromosomes[@]}"
do
  while read -r s1; do
      while read -r s2; do
          [[ "$s1" < "$s2" ]] && echo -e "$s1\t$s2"
      done < /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/contiguous_HORs/HPRC_release2_contiguous_HORs_${chr}.txt
  done < /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/contiguous_HORs/HPRC_release2_contiguous_HORs_${chr}.txt > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/${chr}/HPRC_release2_contiguous_HOR_all_pairs_combinations_${chr}.txt
done
```
Get number of combinations for each chrom
```sh
for chr in "${chromosomes[@]}"
do
  wc -l /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/${chr}/HPRC_release2_contiguous_HOR_all_pairs_combinations_${chr}.txt
done
```
```
80601 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr16/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr16.txt
9591 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr17/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr17.txt
8385 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr18/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr18.txt
76636 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr19/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr19.txt
58653 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr20/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr20.txt
37675 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr21/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr21.txt
55278 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr22/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr22.txt
49455 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chrX/HPRC_release2_contiguous_HOR_all_pairs_combinations_chrX.txt
1176 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chrY/HPRC_release2_contiguous_HOR_all_pairs_combinations_chrY.txt
```
Slurm parameters:
```
#SBATCH --job-name=chr16_pairwise-centrolign
#SBATCH --array=[1-40000]%128
#SBATCH --array=[40001-80601]%128
CHR=chr16

#SBATCH --job-name=chr17_pairwise-centrolign
#SBATCH --array=[1-9591]%128
CHR=chr17

#SBATCH --job-name=chr18_pairwise-centrolign
#SBATCH --array=[1-8385]%128
CHR=chr18

#SBATCH --job-name=chr19_pairwise-centrolign
#SBATCH --array=[1-30000]%128
#SBATCH --array=[30001-76636]%128
CHR=chr19

#SBATCH --job-name=chr20_pairwise-centrolign
#SBATCH --array=[1-30000]%128
#SBATCH --array=[30001-58653]%128
CHR=chr20

#SBATCH --job-name=chr21_pairwise-centrolign
#SBATCH --array=[1-37675]%128
CHR=chr21

#SBATCH --job-name=chr22_pairwise-centrolign
#SBATCH --array=[1-55278]%128
CHR=chr22

#SBATCH --job-name=chrX_pairwise-centrolign
#SBATCH --array=[1-49455]%128
CHR=chrX

#SBATCH --job-name=chrY_pairwise-centrolign
#SBATCH --array=[1-1176]%128
CHR=chrY
```

Slurm script for running pairwise centrolign
```sh
#!/bin/bash
# Slurm script to use centrolign as a pairwise aligner on all pairs of sequences
# from an input directory

#SBATCH --job-name=chrY_pairwise-centrolign
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=56gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=[1-1176]%128
#SBATCH --output=logs/array_job_%A_task_%a.log
#SBATCH --time=1:00:00

date
hostname
pwd

CHR=chrY
CHROMDIR=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/${CHR}
FASTADIR=/private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/
WORKDIR=$CHROMDIR/work/
OUTDIR=$CHROMDIR/pairwise_cigar/

mkdir -p $OUTDIR
mkdir -p $WORKDIR
mkdir -p $FASTADIR

cd $WORKDIR

SAMPLE1=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$CHROMDIR"/HPRC_release2_contiguous_HOR_all_pairs_combinations_${CHR}.txt | cut -f1)
SAMPLE2=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$CHROMDIR"/HPRC_release2_contiguous_HOR_all_pairs_combinations_${CHR}.txt | cut -f2)
echo "sample 1:" $SAMPLE1
echo "sample 2:" $SAMPLE2
echo "out:" $OUTDIR

FASTA1=`grep ${SAMPLE1} /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/HPRC_release2_contiguous_HORs_${CHR}.fasta_list.txt`
FASTA1="/private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/${FASTA1}"
FASTA2=`grep ${SAMPLE2} /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/HPRC_release2_contiguous_HORs_${CHR}.fasta_list.txt`
FASTA2="/private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/${FASTA2}"

echo "fasta 1:" $FASTA1
echo "fasta 2:" $FASTA2
TEMP_FASTA=${WORKDIR}/${SAMPLE1}_${SAMPLE2}.fa
cat $FASTA1 $FASTA2 > $TEMP_FASTA
time /private/home/mmastora/progs/centrolign/build/centrolign -v 3 --skip-calibration $TEMP_FASTA > $OUTDIR/pairwise_cigar_${SAMPLE1}_${SAMPLE2}.txt
rm $TEMP_FASTA
```
#### Test version of locate_hors_from_censat with relaxed discontiguity filters

For each chr, get list of samples with contiguous arrays
```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC_relaxed/release2
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC_relaxed/release2/contiguous_HORs/

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

for chr in "${chromosomes[@]}"
do
    echo "Processing $chr"
    cut -f 1 -d"," extract_hors_HPRC_relaxed_release2.csv | grep -v "sample_id" | \
    while read line ; do
        ls $line/analysis/extract_hors_HPRC_relaxed_outputs/*${chr}_hor_array.fasta | cut -f4 -d"/" | cut -f3 -d"_"
      done | grep -v "cannot access" > /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC_relaxed/release2/contiguous_HORs/HPRC_release2_contiguous_HORs_${chr}.txt
done

# Count number of contiguous arrays per chromosome
for chr in "${chromosomes[@]}"
do  
  count=`wc -l contiguous_HORs/HPRC_release2_contiguous_HORs_${chr}.txt | cut -f1 -d" "`
  echo $chr $count
done
```
```
chr1 398
chr2 442
chr3 195
chr4 221
chr5 375
chr6 222
chr7 369
chr8 373
chr9 424
chr10 385
chr11 433
chr12 401
chr13 303
chr14 348
chr15 358
chr16 417
chr17 158
chr18 150
chr19 392
chr20 354
chr21 275
chr22 333
chrX 324
chrY 60
```

Copy all bed files into one directory:
```
cut -f 1 -d"," extract_hors_HPRC_relaxed_release2.csv | grep -v "sample_id" | \
while read line ; do
  cp ${line}/analysis/extract_hors_HPRC_relaxed_outputs/${line}_hor_arrays.bed contiguous_HORs_bed_files/
done
```
Checking contigs in disagreement with hailey's script
```sh
grep -v "Contig" /private/groups/patenlab/mira/contigs_disagreement.txt | while read line ; do
  SAMPLE=`echo $line | cut -f1 -d"#"`
  CONTIG=`echo $line | cut -f2-3 -d"#"`
  awk -v target="$CONTIG" '
      /^filtering/ { filter_line = $0 }
      {
          if ($0 ~ target) {
              if ($0 ~ /filtering/) {
                  print $0
              } else {
                  print filter_line
                  print $0
              }
          }
      }
  ' /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC_relaxed/release2/${SAMPLE}*/analysis/extract_hors_HPRC_relaxed_outputs/${SAMPLE}*_locate_hors_from_censat.log
done > /private/groups/patenlab/mira/contigs_disagreement_reasons.log

grep "filtering" /private/groups/patenlab/mira/contigs_disagreement_reasons.log | sort | uniq -c
```

### Add CHM13 to release 2 all pairs list

> Adding a .1 to CHM13 label even though it is haploid, for compatibility with downstream scripts

Extract CHM13 alpha arrays
```sh
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

cd /private/groups/patenlab/mira/centrolign/annotations/chm13

mkdir -p per_chrom_for_centrolign/work

for chr in "${chromosomes[@]}"
do
  # bedtools intersect alpha annotation with censat annotation to just get active hors
  echo "processing ${chr} "
  grep -w $chr /private/groups/patenlab/mira/centrolign/annotations/chm13/chm13v2.0_censat_v2.1.bed | grep "hor" | grep -v "dhor" > /private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom_for_centrolign/work/chm13v2.0_censat_v2.1.${chr}.hor.bed

  docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups \
  pegi3s/bedtools bedtools intersect -wa \
      -a /private/groups/patenlab/mira/centrolign/annotations/chm13/chm13v2.0.labels.as_hor.bed \
      -b /private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom_for_centrolign/work/chm13v2.0_censat_v2.1.${chr}.hor.bed \
      > /private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom_for_centrolign/work/chm13v2.0.labels.as_hor.${chr}.active.bed

  # get start and end of active array from as_hor bed
  awk 'NR==1{start=$2; chrom=$1} END{print chrom"\t"start"\t"$3}' \
  /private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom_for_centrolign/work/chm13v2.0.labels.as_hor.${chr}.active.bed \
  > /private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom_for_centrolign/work/chm13v2.0.labels.as_hor.${chr}.active.mrg.bed

  # extract fasta sequence for active hor array
  docker run -it -u `id -u`:`id -g` -v /private/groups:/private/groups \
  pegi3s/bedtools bedtools getfasta \
    -fi /private/groups/patenlab/mira/data/chm13v2.0.fa \
    -bed /private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom_for_centrolign/work/chm13v2.0.labels.as_hor.${chr}.active.mrg.bed | sed '/^>/b; s/[a-z]/\U&/g' | sed '/^>/s/:.*//' | sed "s/>/>CHM13.1 /g" \
    > /private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom_for_centrolign/chm13v2.0.${chr}.active_hor.upper.fa

done
```
Generate lists of CHM13 and every other sample in release 2
```sh
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")


for chr in "${chromosomes[@]}"
do
  while read -r s1; do
      echo -e "CHM13.1\t$s1"
done < /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/contiguous_HORs/HPRC_release2_contiguous_HORs_${chr}.txt > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/${chr}/HPRC_release2_contiguous_HOR_CHM13_combinations_${chr}.txt
done
```

Get number of combinations for each chrom
```sh
for chr in "${chromosomes[@]}"
do
  wc -l /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/${chr}/HPRC_release2_contiguous_HOR_CHM13_combinations_${chr}.txt
done

398 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr1/HPRC_release2_contiguous_HOR_CHM13_combinations_chr1.txt
433 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr2/HPRC_release2_contiguous_HOR_CHM13_combinations_chr2.txt
160 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr3/HPRC_release2_contiguous_HOR_CHM13_combinations_chr3.txt
196 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr4/HPRC_release2_contiguous_HOR_CHM13_combinations_chr4.txt
375 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr5/HPRC_release2_contiguous_HOR_CHM13_combinations_chr5.txt
191 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr6/HPRC_release2_contiguous_HOR_CHM13_combinations_chr6.txt
337 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr7/HPRC_release2_contiguous_HOR_CHM13_combinations_chr7.txt
358 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr8/HPRC_release2_contiguous_HOR_CHM13_combinations_chr8.txt
407 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr9/HPRC_release2_contiguous_HOR_CHM13_combinations_chr9.txt
356 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr10/HPRC_release2_contiguous_HOR_CHM13_combinations_chr10.txt
422 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr11/HPRC_release2_contiguous_HOR_CHM13_combinations_chr11.txt
385 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr12/HPRC_release2_contiguous_HOR_CHM13_combinations_chr12.txt
303 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr13/HPRC_release2_contiguous_HOR_CHM13_combinations_chr13.txt
348 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr14/HPRC_release2_contiguous_HOR_CHM13_combinations_chr14.txt
331 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr15/HPRC_release2_contiguous_HOR_CHM13_combinations_chr15.txt
402 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr16/HPRC_release2_contiguous_HOR_CHM13_combinations_chr16.txt
139 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr17/HPRC_release2_contiguous_HOR_CHM13_combinations_chr17.txt
130 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr18/HPRC_release2_contiguous_HOR_CHM13_combinations_chr18.txt
392 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr19/HPRC_release2_contiguous_HOR_CHM13_combinations_chr19.txt
343 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr20/HPRC_release2_contiguous_HOR_CHM13_combinations_chr20.txt
275 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr21/HPRC_release2_contiguous_HOR_CHM13_combinations_chr21.txt
333 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr22/HPRC_release2_contiguous_HOR_CHM13_combinations_chr22.txt
315 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chrX/HPRC_release2_contiguous_HOR_CHM13_combinations_chrX.txt
49 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chrY/HPRC_release2_contiguous_HOR_CHM13_combinations_chrY.txt
```

Run CHM13 all pairs - adding to pairwise dirs
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

sbatch \
    --job-name=chr1_all_pairs_CHM13 \
    --array=[1-398]%128 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2/scripts/centrolign_all_pairs_CHM13.sh \
    --chr chr1

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr2

sbatch \
    --job-name=chr2_all_pairs_CHM13 \
    --array=[1-433]%128 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2/scripts/centrolign_all_pairs_CHM13.sh \
    --chr chr2

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr3

sbatch \
    --job-name=chr3_all_pairs_CHM13 \
    --array=[1-160]%128 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2/scripts/centrolign_all_pairs_CHM13.sh \
    --chr chr3

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr4

sbatch \
    --job-name=chr4_all_pairs_CHM13 \
    --array=[1-196]%32 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2/scripts/centrolign_all_pairs_CHM13.sh \
    --chr chr4

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr5

sbatch \
    --job-name=chr5_all_pairs_CHM13 \
    --array=[1-375]%32 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2/scripts/centrolign_all_pairs_CHM13.sh \
    --chr chr5
```
