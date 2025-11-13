#!/bin/bash
#SBATCH --job-name=pairwise_consistency
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=7-00:00
#SBATCH --array=[23-26]%5

CSV_FILE=/private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/centrolign_results.csv
CHR_SUBGROUP=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$CSV_FILE" | cut -f1 -d",")
CHR=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$CSV_FILE" | cut -f1 -d"," | cut -f1 -d"_")
INDUCED_CIGARS=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$CSV_FILE" | cut -f10 -d",")

echo $CHR_SUBGROUP

time python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/pairwise_consistency.py \
    ${INDUCED_CIGARS}/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/${CHR}/pairwise_cigar/pairwise_cigar_ \
    > /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/HPRC_r2_QCv2_${CHR_SUBGROUP}_pairwise_consistency.txt
