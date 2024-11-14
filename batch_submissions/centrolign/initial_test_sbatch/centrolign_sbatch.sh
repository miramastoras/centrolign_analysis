#!/bin/bash
#SBATCH --job-name=centrolign_initial
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=500gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm_logs/centrolign_initial_test_sbatch_%x_%j_%A_%a.log
#SBATCH --time=72:00:00
#SBATCH --array=[6-7]%2

set -eux -o pipefail
set -o xtrace

SAMPLE_CSV=$1

# Read the CSV file and extract the sample ID for the current job array task
# Skip first row to avoid the header
SAMPLE_ID=$(awk -F ',' -v task_id=${SLURM_ARRAY_TASK_ID} 'NR>1 && NR==task_id+1 {print $1}' "${SAMPLE_CSV}")
HOR_FASTA=$(awk -F ',' -v task_id=${SLURM_ARRAY_TASK_ID} 'NR>1 && NR==task_id+1 {print $2}' "${SAMPLE_CSV}")
TREE=$(awk -F ',' -v task_id=${SLURM_ARRAY_TASK_ID} 'NR>1 && NR==task_id+1 {print $3}' "${SAMPLE_CSV}")

CHR=$SAMPLE_ID
SAVE_FILES=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/${CHR}/jobstore/
mkdir -p ${SAVE_FILES}

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/${CHR}/

/private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S ${SAVE_FILES} \
    -T ${TREE} \
    ${HOR_FASTA} \
    > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/${CHR}/initial_test_${CHR}.centrolign.gfa
