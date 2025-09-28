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


S3_FILE=/private/groups/patenlab/mira/centrolign/annotations/assemblies_pre_release_v0.6.1.index.csv

S3_ASM=$(grep $ASM_ID $S3_FILE | cut -f13 -d",")
S3_ASM_FAI=$(grep $ASM_ID $S3_FILE | cut -f11 -d",")

LOCAL_FOLDER=/data/tmp/$(whoami)/${CHR}/
mkdir -p ${LOCAL_FOLDER}

# Download assembly
aws s3 cp $S3_ASM $LOCAL_FOLDER/
aws s3 cp $S3_ASM $LOCAL_FOLDER/

ASM_FASTA=$(basename $S3_ASM)

HOR_ARRAY_BED=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/${SMP}.${HAP}_asat_arrays.bed

for CHR in {1..22} X Y M; do
    echo "chr${CHR}"
    OUTPATH=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/chr${CHR}
    mkdir -p $OUTPATH

    #REGIONFILE=/data/tmp/$(whoami)/${CHR}/${ASM_ID}.chr${CHR}.hor.txt
    REGIONFILE=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/test/${ASM_ID}.chr${CHR}.hor.txt
    touch $REGIONFILE
    grep -w chr${CHR} ${HOR_ARRAY_BED} | awk '{ printf "%s:%d-%d\n", $1, $2+1, $3 }' > ${REGIONFILE}
    ls $REGIONFILE

    if [ -s $REGIONFILE ];
        then
            echo "chr${CHR} exists, $REGIONFILE"
            # extract and add the sample name as the sequence name
            echo "extract region" `cat $REGIONFILE`

            HORFASTA=$OUTPATH/${ASM_ID}_${SMP}.${HAP}_chr${CHR}_hor_array.fasta

            samtools faidx -r $REGIONFILE $LOCAL_FOLDER/$ASM_FASTA | sed "s/>/>$SMP.$HAP /g" > $HORFASTA

    else
        echo "~{sampleID} chr${CHR} was filtered out"
    fi
done
