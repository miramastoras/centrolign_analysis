#!/bin/bash
#SBATCH --job-name=annotate_sat_neighbors
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --time=0:30:00
#SBATCH --output=logs/annotate_sat_neighbors_%A_%a.out
#
# For each gene TSV, appends two columns:
#   between_satellites_100kb : flanked on both sides by censat blocks > 100 kb
#   on_acrocentric_short_arm : acrocentric chrom + gene left of active_hor start
#
# Setup (run once before sbatch):
#   bash annotate_satellite_neighbors.sh --make-jobs
#
# Submit:
#   sbatch --array=[1-N]%128 annotate_satellite_neighbors.sh

set -euo pipefail

TSV_DIR='/private/groups/cgl/pnhebbar/pga/censat_gene/intersect_censat_by_chr'
BED_DIR='/private/groups/migalab/juklucas/censat_regions/censat_beds'
OUT_DIR='/private/groups/patenlab/mira/centrolign/analysis/Fig1_genes/satellite_neighbors'
REGIONS_DIR='/private/groups/migalab/juklucas/censat_regions/censat_arrays/tables/censat_regions_pass_qc_by_chrom'
JOBS_TSV="${OUT_DIR}/jobs.tsv"
SCRIPT_DIR="/private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/"

# ── generate jobs list ────────────────────────────────────────────────────────
if [[ "${1:-}" == "--make-jobs" ]]; then
    mkdir -p "$OUT_DIR"
    rm -f "$JOBS_TSV"
    n=0; skip=0
    for tsv in "$TSV_DIR"/*.tsv; do
        fname=$(basename "$tsv")
        sample=$(echo "$fname" | cut -d_ -f1)
        hap=$(echo "$fname"    | cut -d_ -f2)
        bed=$(ls "$BED_DIR"/${sample}_${hap}_*.cenSat.bed 2>/dev/null | head -1 || true)
        if [ -z "$bed" ]; then
            echo "No BED: $sample $hap — $fname" >&2
            (( skip++ )) || true
            continue
        fi
        echo -e "$tsv\t$bed" >> "$JOBS_TSV"
        (( n++ )) || true
    done
    echo "Jobs written: $n  (skipped: $skip)  →  $JOBS_TSV"
    echo "Array range: 1-$n"
    exit 0
fi

# ── per-task processing ───────────────────────────────────────────────────────
source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate base

mkdir -p "$OUT_DIR"

LINE=$(awk -v i="${SLURM_ARRAY_TASK_ID}" 'NR==i' "$JOBS_TSV")
TSV=$(echo "$LINE" | cut -f1)
BED=$(echo "$LINE" | cut -f2)

fname=$(basename "$TSV")
out="${OUT_DIR}/${fname/_intersect.tsv/_sat_neighbors.tsv}"

echo "Task ${SLURM_ARRAY_TASK_ID}: $fname"

python3 "${SCRIPT_DIR}/annotate_satellite_neighbors.py" "$TSV" "$BED" "$out" --regions-dir "$REGIONS_DIR"
