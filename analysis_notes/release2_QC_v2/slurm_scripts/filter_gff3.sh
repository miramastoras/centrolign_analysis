#!/bin/bash
#SBATCH --job-name=filter_gff3_censat
#SBATCH --array=2%463
#SBATCH --cpus-per-task=1
#SBATCH --mem=56gb
#SBATCH --time=1:00:00
#SBATCH --output=logs/filter_gff3_censat_%A_%a.out

set -euo pipefail

source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate base

OUTDIR=$1
JOBS_TSV="${OUTDIR}/jobs.tsv"
RESULTS_DIR="${OUTDIR}/results"
mkdir -p "${RESULTS_DIR}" logs

# Read this job's row (skip header)
LINE=$(awk -v i=$((SLURM_ARRAY_TASK_ID + 2)) 'NR==i' "${JOBS_TSV}")
SAMPLE_ID=$(echo "$LINE" | cut -f1)
HAPLOTYPE=$(echo "$LINE" | cut -f2)
GFF3=$(echo "$LINE" | cut -f3)
BED=$(echo "$LINE" | cut -f4)

OUT="${RESULTS_DIR}/${SAMPLE_ID}.${HAPLOTYPE}.genes.gff3"

echo "Processing ${SAMPLE_ID} hap${HAPLOTYPE}"

# Subset GFF3 to gene features overlapping censat regions,
# then drop any gene whose name starts with ENS
awk 'BEGIN{OFS="\t"} {split($1,a,"#"); $1=a[3]; print}' "${BED}" \
    | bedtools intersect \
        -a "${GFF3}" \
        -b stdin \
        -wo \
    | awk '$3 == "gene"' \
    | awk '
        {
            chrom = $13
            overlap = $NF
            name = ""
            n = split($9, attrs, ";")
            for (i = 1; i <= n; i++) {
                if (attrs[i] ~ /^source_gene_common_name=/) {
                    split(attrs[i], kv, "="); name = kv[2]; break
                }
                if (attrs[i] ~ /^gene_name=/ && name == "") {
                    split(attrs[i], kv, "="); name = kv[2]
                }
                if (attrs[i] ~ /^Name=/ && name == "") {
                    split(attrs[i], kv, "="); name = kv[2]
                }
            }
            if (name !~ /^ENS/)
                print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"chrom"\t"overlap
        }
    ' \
    > "${OUT}"

echo "Done: $(wc -l < ${OUT}) genes written to ${OUT}"
