#!/bin/bash
# Slurm script to run the sim_centromere script to generate pairwise sequence
# alignment problems
#SBATCH --job-name=make-pair-sim-cases
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=100gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-60
#SBATCH --output=make_sim_cases_slurm_logs/array_job_%A_task_%a.log
#SBATCH --time=12:00:00

date
hostname
pwd
echo "array job" $SLURM_ARRAY_TASK_ID

# note: sample size is handled in sbatch array size
CHR=$1
DATE=20250421

SIMDIR=/private/groups/patenlab/mira/centrolign/simulations/pairwise_simulations/
OUTPARDIR=$SIMDIR/pair_"$CHR"_sim_cases_"$DATE"/

SIM_CENTROMERE=/private/home/mmastora/progs/centrolign/build/sim_centromere

# the base array that we'll simulate
FASTA=/private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom/chm13v2.0.${CHR}.active_hor.upper.fa
BED=/private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom/chm13v2.0.labels.as_hor.${CHR}.active.shifted.bed

GENS=("25" "50" "100" "150" "200" "300")
i=$(($SLURM_ARRAY_TASK_ID % 6))
GEN=${GENS[$i]}
OUTDIR=$OUTPARDIR/gen"$GEN"/case_"$SLURM_ARRAY_TASK_ID"/
mkdir -p $OUTDIR
# simulate the sequence
echo "generating sequence pair" $SLURM_ARRAY_TASK_ID "with" $GEN "generations"
/usr/bin/time -v $SIM_CENTROMERE -g $GEN -o $OUTDIR/sim $FASTA $BED
# record this case for the centrolign file
echo gen"$GEN"/case_"$SLURM_ARRAY_TASK_ID" >> $OUTPARDIR/cases.txt

# this is silly, but it makes the evaluation logic much simpler if we always have a tree
printf "(seq1,seq2):%d;\n" $GEN > $OUTDIR/tree.txt
cat $OUTDIR/sim*.fasta > $OUTDIR/all_seqs.fasta
