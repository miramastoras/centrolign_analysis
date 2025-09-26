#!/bin/bash
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

#510
source /private/groups/patenlab/jeizenga/centromere/venv/bin/activate

COMBINATIONS_FILE=/private/groups/patenlab/mira/centrolign/simulations/centrolign_pairwise_vs_MSA/combination_lists/all_chroms_all_cases.txt

CHR=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$COMBINATIONS_FILE" | cut -f1 -d",")
CASE=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$COMBINATIONS_FILE" | cut -f2 -d",")

echo $CHR, $CASE
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/analyze_case_v2.py \
  /private/groups/patenlab/mira/centrolign/simulations/MSA_simulations/msa_${CHR}_sim_cases_20250402/${CASE}/ \
  /private/home/mmastora/progs/centrolign/build/tree_pair_dist \
  /private/home/mmastora/progs/centrolign/build/compare_truth_aln \
  /private/groups/patenlab/mira/centrolign/simulations/centrolign_pairwise_vs_MSA/pairwise_cigars/${CHR}/${CASE}/pairwise_cigar/ \
  /private/groups/patenlab/mira/centrolign/simulations/centrolign_pairwise_vs_MSA/analyze_case_results/
