#!/bin/bash
#SBATCH --job-name=extract_hor_fastas
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=56gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=[1-127]%127
#SBATCH --exclude=phoenix-[09,10,22,23,24,18]
#SBATCH --output=logs/array_job_%A_task_%a.log
#SBATCH --time=1:00:00

SMP_FILE=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_samples_with_asats.txt

SMP=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$SMP_FILE" | cut -f1 -d",")
HAP=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$SMP_FILE" | cut -f2 -d",")
ASM_ID=${SMP}_${HAP}
HAP_CODE=$(echo "$HAP" | awk '{print ($1 == "hap1") ? 1 : 2}')

S3_FILE=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/HGSVC_assemblies.csv

S3_ASM=$(grep $ASM_ID $S3_FILE | cut -f2 -d",")
S3_ASM_FAI=$(grep $ASM_ID $S3_FILE | cut -f3 -d",")

LOCAL_FOLDER=/data/tmp/$(whoami)/${CHR}/
mkdir -p ${LOCAL_FOLDER}

# Download assembly
wget -P $LOCAL_FOLDER/ $S3_ASM 
wget -P $LOCAL_FOLDER/ $S3_ASM_FAI

ASM_FASTA=$(basename $S3_ASM)

HOR_ARRAY_BED=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_per_smp_asat_beds/${SMP}.${HAP}_asat_arrays.bed

for CHR in {1..22} X Y; do
    echo "chr${CHR}"
    OUTPATH=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/extract_fastas_hgsvc/chr${CHR}
    mkdir -p $OUTPATH

    #REGIONFILE=/data/tmp/$(whoami)/${CHR}/${ASM_ID}.chr${CHR}.hor.txt
    REGIONFILE=${LOCAL_FOLDER}/${ASM_ID}.chr${CHR}.hor.txt
    touch $REGIONFILE
    grep -w chr${CHR} ${HOR_ARRAY_BED} | awk '{ printf "%s:%d-%d\n", $1, $2+1, $3 }' > ${REGIONFILE}
    ls $REGIONFILE

    if [ -s $REGIONFILE ];
        then
            echo "chr${CHR} exists, $REGIONFILE"
            # extract and add the sample name as the sequence name
            echo "extract region" `cat $REGIONFILE`

            HORFASTA=$OUTPATH/${SMP}.${HAP_CODE}_chr${CHR}_hor_array.fasta

            samtools faidx -r $REGIONFILE $LOCAL_FOLDER/$ASM_FASTA | sed "s/>/>$SMP.$HAP_CODE /g" > $HORFASTA

    else
        echo "chr${CHR} was filtered out"
    fi
done

rm $LOCAL_FOLDER/${ASM_FASTA}
rm $LOCAL_FOLDER/${ASM_FASTA}.fai
