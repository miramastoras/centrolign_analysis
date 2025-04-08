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

```

Slurm script for running pairwise centrolign
```sh
#!/bin/bash
# Slurm script to use centrolign as a pairwise aligner on all pairs of sequences
# from an input directory

#SBATCH --job-name=chr19_pairwise-centrolign
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=56gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=[1-30000]%128
#SBATCH --output=logs/array_job_%A_task_%a.log
#SBATCH --time=1:00:00

date
hostname
pwd

CHR=chr19
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
#### Preparing inputs for centrolign giraffe testing

Identify samples that align well to other samples in Faith's test set GFA, but are not in the GFA

```sh
# all pairwise alignments in release 2 that align well
awk -F, '$3 < 0.5' /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr12/pairwise_distance.csv  > /private/groups/patenlab/mira/centrolign/giraffe/pairwise_dist_lt0.5.csv

# all pairwise alignments involving samples inside Faith's test graph
cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/fasta_list.all_sample_ids.txt | while read line ; do grep $line /private/groups/patenlab/mira/centrolign/giraffe/pairwise_dist_lt0.5.csv ; done | sort | uniq > /private/groups/patenlab/mira/centrolign/giraffe/pairwise_alns_inside_graph.csv

# get all sample names involved in these pairwise alignments
cut -f 1 -d"," /private/groups/patenlab/mira/centrolign/giraffe/pairwise_alns_inside_graph.csv > samps
cut -f 2 -d"," /private/groups/patenlab/mira/centrolign/giraffe/pairwise_alns_inside_graph.csv >> samps

sort samps | uniq > tmp ; mv tmp samps

# now sample names that are NOT in graph
grep -v -f /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/fasta_list.all_sample_ids.txt samps > /private/groups/patenlab/mira/centrolign/giraffe/samples_aligning_well_not_in_graph.txt

# get their fasta locations
cat /private/groups/patenlab/mira/centrolign/giraffe/samples_aligning_well_not_in_graph.txt | while read line ; do grep $line /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/HPRC_release2_contiguous_HORs_chr12.fasta_list.txt ; done > /private/groups/patenlab/mira/centrolign/giraffe/fastas_not_in_graph_to_align.txt
```
