#!/bin/bash
# Slurm script to use centrolign as a pairwise aligner on all pairs of sequences
# from an input directory

#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=56gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/array_job_%A_task_%a.log
#SBATCH --time=1:00:00

date
hostname
pwd

COMBINATIONS_FILE=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/${CHR}/release2_QC_v2_all_pairs_combinations_${CHR}.txt

CHROMDIR=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/${CHR}

WORKDIR=$CHROMDIR/work/
OUTDIR=$CHROMDIR/pairwise_cigar/

mkdir -p $OUTDIR
mkdir -p $WORKDIR

cd $WORKDIR

SAMPLE1=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$COMBINATIONS_FILE" | cut -f1 | xargs basename | cut -f1 -d"_")
SAMPLE2=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$COMBINATIONS_FILE" | cut -f2 | xargs basename | cut -f1 -d"_")
echo "sample 1:" $SAMPLE1
echo "sample 2:" $SAMPLE2
echo "out:" $OUTDIR

FASTA1=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$COMBINATIONS_FILE" | cut -f1)

FASTA2=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$COMBINATIONS_FILE" | cut -f2)

echo "fasta 1:" $FASTA1
echo "fasta 2:" $FASTA2

TEMP_FASTA=${WORKDIR}/${SAMPLE1}_${SAMPLE2}.fa
cat $FASTA1 $FASTA2 > $TEMP_FASTA
time /private/home/mmastora/progs/centrolign/build/centrolign -v 3 --skip-calibration $TEMP_FASTA > $OUTDIR/pairwise_cigar_${SAMPLE1}_${SAMPLE2}.txt
rm $TEMP_FASTA
