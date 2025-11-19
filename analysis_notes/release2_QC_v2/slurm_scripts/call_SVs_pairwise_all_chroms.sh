#!/bin/bash
#SBATCH --job-name=SVs_all_chroms
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/call_SVs_%x.%j.log
#SBATCH --time=7-00:00
#SBATCH --array=[1-28]%28

CLADE_NAMES=$1
OUTDIR=$2
CLADE=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$CLADE_NAMES" | cut -f1 -d",")
CIGAR_PATH=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$CLADE_NAMES" | cut -f2 -d",")
SAMPLE_LIST=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$CLADE_NAMES" | cut -f3 -d",")
CHR=$(echo $CLADE | cut -f1 -d"_")

echo $CHR
echo $CLADE
echo $CIGAR_PATH
echo $SAMPLE_LIST

mkdir -p ${OUTDIR}/${CHR}/SV_beds/${CLADE}/

time python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/call_SVs_pairwise.py \
  -c ${CIGAR_PATH}/pairwise_cigar_ \
  -s ${SAMPLE_LIST} \
  -o ${OUTDIR}/${CHR}/SV_beds/${CLADE}/
