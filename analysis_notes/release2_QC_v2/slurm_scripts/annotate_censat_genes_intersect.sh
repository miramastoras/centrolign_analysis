#!/bin/bash
#SBATCH --job-name=annotate_censat
#SBATCH --array=1-4072%128
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --time=0:30:00
#SBATCH --output=logs/annotate_censat_%A_%a.out

# Annotate per-chrom gene TSVs with overlapping censat categories and % overlap.
# Adds a final column: "active_hor(82.3%);CT(17.7%)" or "." if no overlap.
#
# Setup (run once before sbatch):
#   bash annotate_censat_genes_intersect.sh --make-jobs
#
# Submit:
#   sbatch annotate_censat_genes_intersect.sh

set -euo pipefail

TSV_DIR='/private/groups/cgl/pnhebbar/pga/censat_gene/intersect_censat_by_chr'
BED_DIR='/private/groups/migalab/juklucas/censat_regions/censat_beds'
OUT_DIR='/private/groups/patenlab/mira/centrolign/analysis/Fig1_genes/annotate_censat_prajnas_list'
JOBS_TSV="${OUT_DIR}/jobs.tsv"

# ── generate jobs list (run with --make-jobs before submitting) ───────────────
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
echo "Task ${SLURM_ARRAY_TASK_ID}: $fname"

TMP=$(mktemp -d)
trap "rm -rf $TMP" EXIT

# strip "sample#hap#" prefix so contig names match the plain accessions in the TSV
grep -v '^track' "$BED" \
    | awk 'BEGIN{OFS="\t"} { n=split($1,a,"#"); $1=(n==3 ? a[3] : $1); print }' \
    > "$TMP/censat.bed"

# bedtools intersect -wao:
#   A = TSV (11 cols), B = censat BED (9 cols), last col = overlap bp
#   BED annot (col 4) lands at field 11+4 = 15; overlap = field 21 ($NF)
# awk aggregates multiple overlaps per gene into one semicolon-separated string
bedtools intersect -a "$TSV" -b "$TMP/censat.bed" -wao \
| awk 'BEGIN { OFS="\t" } {
    key = $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11
    ov       = int($NF)
    gene_len = int($5) - int($4)
    if (ov > 0 && gene_len > 0) {
        pct   = 100 * ov / gene_len
        entry = $15 "(" sprintf("%.1f", pct) "%)"
        if (data[key] == "") data[key] = entry
        else                 data[key] = data[key] ";" entry
    } else if (!(key in data)) {
        data[key] = "."
    }
    if (!(key in order)) { order[key] = ++n; keys[n] = key }
}
END {
    for (i = 1; i <= n; i++) print keys[i] "\t" data[keys[i]]
}' \
> "$OUT_DIR/${fname/_intersect.tsv/_censat_annot.tsv}"

echo "Done: $(wc -l < "$OUT_DIR/${fname/_intersect.tsv/_censat_annot.tsv}") records"
