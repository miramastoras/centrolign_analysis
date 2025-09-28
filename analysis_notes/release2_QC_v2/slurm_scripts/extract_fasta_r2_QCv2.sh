#!/bin/bash
#SBATCH --job-name=extract_hor_fastas
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=56gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=[1-2]%2
#SBATCH --output=logs/array_job_%A_task_%a.log
#SBATCH --time=1:00:00

SMP_FILE=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_samples_with_asats.txt

SMP=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$SMP_FILE" | cut -f1 -d",")
HAP=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$SMP_FILE" | cut -f2 -d",")
ASM_ID=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$SMP_FILE" | cut -f3 -d",")

echo $SMP, $HAP, $ASM_ID

S3_FILE=/private/groups/patenlab/mira/centrolign/annotations/assemblies_pre_release_v0.6.1.index.csv

S3_ASM=$(grep $ASM_ID $S3_FILE | cut -f13 -d",")
S3_ASM_FAI=$(grep $ASM_ID $S3_FILE | cut -f11 -d",")

LOCAL_FOLDER=/data/tmp/$(whoami)/${CHR}/

echo $S3_ASM, $S3_ASM_FAI
#mkdir -p ${LOCAL_FOLDER}

# Download assembly
#aws s3 cp $S3_ASM $TMP_DIR/
#aws s3 cp $S3_ASM $TMP_DIR/

#samtools faidx -r $REGIONFILE ~{assemblyFasta} | sed "s/>/>$SAMPLE.$HAP /g" > $HORFASTA
