#!/bin/bash
#SBATCH --job-name=extract_hors_initial_test
#SBATCH --cpus-per-task=8
#SBATCH --threads-per-core=1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=16gb
#SBATCH --time=1:00:00
#SBATCH --partition=short
#SBATCH --output=slurm_logs/submission_%x_%j_%A_%a.log
#SBATCH --array=[1-189]%40

set -e
set -x

###############################################################################
##                             Prepare For Run                               ##
###############################################################################

OUTDIR=/private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_from_assemblies_sbatch
SAMPLE_CSV=$1

# Read the CSV file and extract the sample ID for the current job array task
# Skip first row to avoid the header
SAMPLE_ID=$(awk -F ',' -v task_id=${SLURM_ARRAY_TASK_ID} 'NR>1 && NR==task_id+1 {print $1}' "${SAMPLE_CSV}")

# Ensure a sample ID is obtained
if [ -z "${SAMPLE_ID}" ]; then
    echo "Error: Failed to retrieve a valid sample ID. Exiting."
    exit 1
fi

mkdir -p ${OUTDIR}/${SAMPLE_ID}

# get assembly and paf file
ASSEMBLY=$(awk -F ',' -v task_id=${SLURM_ARRAY_TASK_ID} 'NR>1 && NR==task_id+1 {print $2}' "${SAMPLE_CSV}")
PAF=$(awk -F ',' -v task_id=${SLURM_ARRAY_TASK_ID} 'NR>1 && NR==task_id+1 {print $3}' "${SAMPLE_CSV}")
AS_HOR_SF=$(awk -F ',' -v task_id=${SLURM_ARRAY_TASK_ID} 'NR>1 && NR==task_id+1 {print $4}' "${SAMPLE_CSV}")
CENSAT=$(awk -F ',' -v task_id=${SLURM_ARRAY_TASK_ID} 'NR>1 && NR==task_id+1 {print $5}' "${SAMPLE_CSV}")

# Download assembly
mkdir -p ${OUTDIR}/${SAMPLE_ID}/tmp_files/
aws s3 cp ${ASSEMBLY} ${OUTDIR}/${SAMPLE_ID}/tmp_files/

# assembly needs to be unzipped for script
ASSEMBLY_BASENAME=`basename ${ASSEMBLY}`
gunzip ${OUTDIR}/${SAMPLE_ID}/tmp_files/${ASSEMBLY_BASENAME}

ASM_LOCAL="${OUTDIR}/${SAMPLE_ID}/tmp_files/${ASSEMBLY_BASENAME%.gz}"

ls ${CENSAT}
ls ${AS_HOR_SF}
ls ${PAF}
ls ${ASM_LOCAL}

samtools faidx ${ASM_LOCAL}

# run locate hors script
python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/locate_hors_from_censat.py \
    -c ${CENSAT} \
    -a ${AS_HOR_SF} \
    -p ${PAF} \
    -f ${ASM_LOCAL} \
    > ${OUTDIR}/${SAMPLE_ID}/${SAMPLE_ID}_hor_arrays.bed

### for each assembly, split fasta by chromosome

# switch haplotypes labels
SAMPLE==$(echo "$SAMPLE_ID" | sed 's/_hap[12]//')

if [[ "${SAMPLE_ID}" == *"hap1"* ]]; then
      PARNUM=2
      HPRC_PARENT=hap1
else
  PARNUM=1
  HPRC_PARENT=hap2
fi

CHR=chr12

REGIONFILE=${OUTDIR}/${SAMPLE_ID}/${SAMPLE_ID}.chr12.hor.txt

grep $CHR ${OUTDIR}/${SAMPLE_ID}/${SAMPLE_ID}_hor_arrays.bed | awk '{ printf "%s:%d-%d\n", $1, $2+1, $3 }' > "$REGIONFILE"

STRAND=$(grep $CHR ${OUTDIR}/${SAMPLE_ID}/${SAMPLE_ID}_hor_arrays.bed | cut -f 6)

    if [ -s $REGIONFILE ];
    then

        # extract and add the sample name as the sequence name
        echo "extract region" `cat $REGIONFILE`
        echo "strand" $STRAND

        HORFASTA=${OUTDIR}/${SAMPLE_ID}/${SAMPLE_ID}_parnum_${PARNUM}_hor_array.fasta

        samtools faidx -r $REGIONFILE ${ASM_LOCAL} | sed "s/>/>$SAMPLE.$PARNUM /g" > $HORFASTA

        if [ $STRAND = "-" ];
        then
            echo "reverse complementing sequence"
            TMPFILE=$(mktemp)
            /private/groups/patenlab/mira/centrolign/github/centromere-scripts/data_processing_utils/fasta_to_rev_comp.py $HORFASTA > $TMPFILE
            mv $TMPFILE $HORFASTA
        fi
    else
      echo "$SAMPLE_ID $CHR was filtered out"
    fi

rm -rf ${OUTDIR}/${SAMPLE_ID}/tmp_files/
