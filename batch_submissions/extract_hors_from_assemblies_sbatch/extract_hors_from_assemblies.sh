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
#SBATCH --array=[1]%1

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

# run locate hors script
python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/locate_hors_from_censat.py \
    -c ${CENSAT} \
    -a ${AS_HOR_SF} \
    -p ${PAF} \
    -f ${ASM_LOCAL} \
    > ${OUTDIR}/${SAMPLE_ID}/${SAMPLE_ID}_hor_arrays.bed

# for each assembly, split fasta by chromosome
if [[ "${SAMPLE_ID}" == *"$hap1"* ]]; then
  then
      PARNUM=2
else
  PARNUM=1
fi

REGIONFILE=${OUTDIR}/${SAMPLE_ID}/${SAMPLE_ID}_hor_arrays.bed
