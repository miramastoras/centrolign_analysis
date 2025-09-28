#!/bin/bash
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=56gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/array_job_%A_task_%a.log
#SBATCH --time=1:00:00

ALL_SAMPLES=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/active_array_csvs/asat_arrays_${CHR}.csv

SAMPLE=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$QC_FILE" | cut -f1 -d",")
HAP=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$QC_FILE" | cut -f2 -d",")

S3_ASM=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$QC_FILE" | cut -f7 -d",")
S3_ASM_FAI=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$QC_FILE" | cut -f8 -d",")

LOCAL_FOLDER=/data/tmp/$(whoami)/${CHR}/
mkdir -p ${LOCAL_FOLDER}

# Download assembly
aws s3 cp $S3_ASM $TMP_DIR/
aws s3 cp $S3_ASM $TMP_DIR/

samtools faidx -r $REGIONFILE ~{assemblyFasta} | sed "s/>/>$SAMPLE.$HAP /g" > $HORFASTA
