#!/bin/bash
#SBATCH --job-name=annotate_censat_genes
#SBATCH --array=2-463%128
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --time=1:00:00
#SBATCH --output=logs/annotate_censat_%A_%a.out

set -euo pipefail

source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate base

OUTDIR=$1
JOBS_TSV="${OUTDIR}/jobs.tsv"
RESULTS_DIR="${OUTDIR}/results"

LINE=$(awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i' "${JOBS_TSV}")
SAMPLE_ID=$(echo "$LINE" | cut -f1)
HAPLOTYPE=$(echo "$LINE" | cut -f2)
GFF3=$(echo "$LINE" | cut -f3)
CENSAT=$(echo "$LINE" | cut -f4)

OUT="${RESULTS_DIR}/${SAMPLE_ID}.${HAPLOTYPE}.genes.censat.tsv"

echo "Processing ${SAMPLE_ID} hap${HAPLOTYPE}"

grep -v "^track" "${CENSAT}" \
    | awk 'BEGIN{OFS="\t"} {split($1,a,"#"); $1=a[3]; print}' \
    | bedtools intersect \
        -a "${GFF3}" \
        -b stdin \
        -wo \
    | awk '$3 == "gene"' \
    | awk 'BEGIN{OFS="\t"} {
        annot   = $14
        overlap = $20

        if (annot ~ /^active_hor/) cat = "active_hor"
        else if (annot ~ /^dhor/) cat = "dhor"
        else if (annot ~ /^hor/) cat = "hor"
        else if (annot ~ /^mixedAlpha/) cat = "mixedAlpha"
        else if (annot ~ /^HSat1|HSAT1/) cat = "HSAT1"
        else if (annot ~ /^HSat2|HSAT2/) cat = "HSAT2"
        else if (annot ~ /^HSat3|HSAT3/) cat = "HSAT3"
        else if (annot ~ /HSAT4/) cat = "HSAT4"
        else if (annot ~ /HSAT5/) cat = "HSAT5"
        else if (annot ~ /^gSat/) cat = "gSat"
        else if (annot ~ /^bSat/) cat = "bSat"
        else if (annot ~ /^mon/) cat = "mon"
        else if (annot ~ /^[Cc][Tt]/) cat = "CT"
        else if (annot ~ /^GAP/) cat = "GAP"
        else if (annot ~ /^rDNA|5SRNA/) cat = "rDNA"
        else if (annot ~ /SATR/) cat = "SATR"
        else if (annot ~ /SST1/) cat = "SST1"
        else if (annot ~ /CER/) cat = "CER"
        else if (annot ~ /ACRO/) cat = "acrocentric"
        else cat = "other"

        print $1,$4,$5,$8,$10,$9,annot,cat,overlap
    }' \
    > "${OUT}"

echo "Done: $(wc -l < ${OUT}) records written to ${OUT}"