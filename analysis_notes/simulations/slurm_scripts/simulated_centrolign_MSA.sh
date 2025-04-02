#!/bin/bash
# Slurm script to benchmark the output of a centrolign multiple sequence alignment
# of simulated sequences, created by the sim_centromere script.
# Calls the analyze_case.py and infer_tree.py scripts.

#SBATCH --job-name=simulated-centrolign-msa
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=NONE
#SBATCH --nodes=1
#SBATCH --mem=160gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-30
#SBATCH --output=simulated_centrolign_slurm_logs/analysis_array_job_%A_task_%a.log
#SBATCH --time=12:00:00

date
hostname
pwd

CHR=$1
DATE=20250402

SIMDIR=/private/groups/patenlab/mira/centrolign/simulations/MSA_simulations
CASEDIR=$SIMDIR/msa_"$CHR"_sim_cases_"$DATE"/

CASE=$(ls $CASEDIR | sort | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

source /private/groups/patenlab/jeizenga/centromere/venv/bin/activate
cd $CASEDIR/$CASE
mkdir -p induced
mkdir -p subprobs

CENTROLIGNDIR=/private/home/mmastora/progs/centrolign/build/
CENTROLIGN=$CENTROLIGNDIR/centrolign
TRUTH_COMPARE=$CENTROLIGNDIR/compare_truth_aln
TREE_COMPARE=$CENTROLIGNDIR/tree_compare
TREE_DIST=$CENTROLIGNDIR/tree_pair_dist
ANALYZE_CASE=/private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/analyze_case.py
INFER_TREE=/private/groups/patenlab/mira/centrolign/github/centromere-scripts/data_exploration/infer_tree.py

#echo "beginning alignment of case" $CASEDIR/$CASE
/usr/bin/time -v $CENTROLIGN -v 3 -T tree.txt -A induced/aln -S subprobs/sp all_seqs.fasta > msa.gfa 2> >( tee err.txt >&2 )

echo "alignment completed, analyzing results"

source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate skbio

time $INFER_TREE induced 1 > inferred_tree.txt
time $TREE_COMPARE tree.txt inferred_tree.txt > tree_comparison.tsv

# the output makes more sense if we give it a non-trivial directory
cd $CASEDIR
time $ANALYZE_CASE $CASE $TREE_DIST $TRUTH_COMPARE
