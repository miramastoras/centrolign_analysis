#!/bin/bash
# This slurm script creates all pairs alignments of the simulated sequences used
# in benchmarking the centrolign MSA (created by make_sim_cases_MSA).
# We do this so we can compare the performance of pairwise alignments induced from the MSA,
# and direct pairwise alignments on the same sequences
#!/bin/bash
# Slurm script to use centrolign as a pairwise aligner on all pairs of sequences
# from an input directory

#SBATCH --job-name=msa_simulations_direct_pairwise
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=56gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=[1-2]%128
#SBATCH --output=logs/array_job_%A_task_%a.log
#SBATCH --time=1:00:00

date
hostname
pwd

COMBINATIONS_FILE=/private/groups/patenlab/mira/centrolign/simulations/centrolign_pairwise_vs_MSA/combination_lists/all_combinations.txt

CHR=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$COMBINATIONS_FILE" | cut -f3)
CASE=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$COMBINATIONS_FILE" | cut -f4)

CHROMDIR=/private/groups/patenlab/mira/centrolign/simulations/centrolign_pairwise_vs_MSA/pairwise_cigars/${CHR}/${CASE}

WORKDIR=$CHROMDIR/work/
OUTDIR=$CHROMDIR/pairwise_cigar/

mkdir -p $OUTDIR
mkdir -p $WORKDIR

cd $WORKDIR

SAMPLE1=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$COMBINATIONS_FILE" | cut -f1 | xargs basename | cut -f2 -d"_")
SAMPLE2=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$COMBINATIONS_FILE" | cut -f2 | xargs basename | cut -f2 -d"_")
echo "sample 1:" $SAMPLE1
echo "sample 2:" $SAMPLE2
echo "out:" $OUTDIR

FASTA1=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$COMBINATIONS_FILE" | cut -f1)

FASTA2=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$COMBINATIONS_FILE" | cut -f2)

echo "fasta 1:" $FASTA1
echo "fasta 2:" $FASTA2

TEMP_FASTA=${WORKDIR}/sim_${SAMPLE1}_${SAMPLE2}.fa
cat $FASTA1 $FASTA2 > $TEMP_FASTA
time /private/home/mmastora/progs/centrolign/build/centrolign -v 3 --skip-calibration $TEMP_FASTA > $OUTDIR/aln_${SAMPLE1}_${SAMPLE2}.txt
rm $TEMP_FASTA
