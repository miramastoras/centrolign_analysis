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
    ls | grep hap | \
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
78913 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr1/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr1.txt
18118 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr6/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr6.txt
63823 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr8/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr8.txt
63118 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr10/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr10.txt
73829 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr12/HPRC_release2_contiguous_HOR_all_pairs_combinations_chr12.txt
```
Slurm parameters:
```
#SBATCH --job-name=chr6_pairwise-centrolign
#SBATCH --array=[1-18118]%32
CHR=chr6

#SBATCH --job-name=chr1_pairwise-centrolign
#SBATCH --array=[1-78913]%128
CHR=chr1

#SBATCH --job-name=chr12_pairwise-centrolign
#SBATCH --array=[1-10000]%128
#SBATCH --array=[10001-73829]%128
CHR=chr12

#SBATCH --job-name=chr10_pairwise-centrolign
#SBATCH --array=[1-63118]%128
CHR=chr10

#SBATCH --job-name=chr8_pairwise-centrolign
#SBATCH --array=[1-63823]%32
CHR=chr8
```
Slurm script for running pairwise centrolign
```sh
#!/bin/bash
# Slurm script to use centrolign as a pairwise aligner on all pairs of sequences
# from an input directory

#SBATCH --job-name=chr8_pairwise-centrolign
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=56gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=[1-12720]%32
#SBATCH --output=logs/array_job_%A_task_%a.log
#SBATCH --time=1:00:00

date
hostname
pwd

CHR=chr8
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

S1_hap=`echo $SAMPLE1 | cut -f2 -d"."`
S2_hap=`echo $SAMPLE2 | cut -f2 -d"."`

S1_ID=`echo $SAMPLE1 | cut -f1 -d"."`
S2_ID=`echo $SAMPLE2 | cut -f1 -d"."`

if [ ${S1_hap} == 1 ];
then
    S1_HPRC_hap="hap1"
else
    S1_HPRC_hap="hap2"
fi

if [ ${S2_hap} == 1 ];
then
    S2_HPRC_hap="hap1"
else
    S2_HPRC_hap="hap2"
fi

FASTA1=${FASTADIR}${S1_ID}_${S1_HPRC_hap}/analysis/extract_hors_HPRC_outputs/${S1_ID}_${S1_HPRC_hap}_${SAMPLE1}_${CHR}_hor_array.fasta

FASTA2=${FASTADIR}${S2_ID}_${S2_HPRC_hap}/analysis/extract_hors_HPRC_outputs/${S2_ID}_${S2_HPRC_hap}_${SAMPLE2}_${CHR}_hor_array.fasta

echo "fasta 1:" $FASTA1
echo "fasta 2:" $FASTA2
TEMP_FASTA=${WORKDIR}/${SAMPLE1}_${SAMPLE2}.fa
cat $FASTA1 $FASTA2 > $TEMP_FASTA
time /private/home/mmastora/progs/centrolign/build/centrolign -v 3 --skip-calibration $TEMP_FASTA > $OUTDIR/pairwise_cigar_${SAMPLE1}_${SAMPLE2}.txt
rm $TEMP_FASTA
```
