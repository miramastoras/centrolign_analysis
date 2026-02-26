### Comparing the mutation rate on the flanks with inside the array

#### Extracting aligned bases from the graph

Get CHM13 +/- 50kb of the flanks
```sh
ls /private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom/work/ | grep .active.mrg.bed | while read line ; do
    cat $line >>/private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom/work/chm13v2.0.labels.as_hor.all_chrs.active.bed
  done

cut -f1,2 /private/groups/patenlab/mira/data/chm13v2.0.fa.fai > /private/groups/patenlab/mira/data/chm13v2.0.chrom_sizes.txt

bedtools flank -i /private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom/work/chm13v2.0.labels.as_hor.all_chrs.active.bed -g /private/groups/patenlab/mira/data/chm13v2.0.chrom_sizes.txt -b 50000 > /private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom/work/chm13v2.0.labels.as_hor.all_chrs.50kb_flanks.bed
```

Convert HPRC HAL file to MAF
```sh
source /private/groups/cgl/cactus/venv-cactus-latest/bin/activate

# concatenate bed regions from CHM13
hal2maf \
  --refGenome CHM13 \
  --refTargets /private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom/work/chm13v2.0.labels.as_hor.all_chrs.50kb_flanks.bed \
  --noAncestors \
  /private/groups/patenlab/mira/centrolign/analysis/mutation_rate_comparison_to_flanks/hprc-v2.0-mc-chm13.full.hal \
  /private/groups/patenlab/mira/centrolign/analysis/mutation_rate_comparison_to_flanks/hprc-v2.0-mc-chm13.50kb_flanks.maf

halStats --genomes  /private/groups/patenlab/mira/centrolign/analysis/mutation_rate_comparison_to_flanks/hprc-v2.0-mc-chm13.full.hal
```
Convert maf to bed file per sample
```sh
python /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/maf_to_bed_per_sample.py \
    --maf /private/groups/patenlab/mira/centrolign/analysis/mutation_rate_comparison_to_flanks/hprc-v2.0-mc-chm13.50kb_flanks.maf \
    --ref CHM13 \
    --outdir /private/groups/patenlab/mira/centrolign/analysis/mutation_rate_comparison_to_flanks/per_sample_maf_beds
```

Sanity check: confirm that missing genotypes in the vcf correspond to places uncovered by the bed files:
```sh
python /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/check_missing_concordance_hal_bed.py \
    --vcf /private/groups/patenlab/mira/centrolign/analysis/mutation_rate_comparison_to_flanks/flanking_50kb/chr20_flanking_50kb.vcf.gz \
    --beddir /private/groups/patenlab/mira/centrolign/analysis/mutation_rate_comparison_to_flanks/per_sample_maf_beds/ \
    --out chr20_concordance.tsv
```

Merge bedfiles to cover cases where CHM13 is traversed multiple times
```sh
cd /private/groups/patenlab/mira/centrolign/analysis/mutation_rate_comparison_to_flanks/per_sample_maf_beds

for bed in /private/groups/patenlab/mira/centrolign/analysis/mutation_rate_comparison_to_flanks/per_sample_maf_beds/*.bed; do
    sample=$(basename "$bed" .bed)
    bedtools merge -i "$bed" > "/private/groups/patenlab/mira/centrolign/analysis/mutation_rate_comparison_to_flanks/per_sample_maf_beds/${sample}.mrg.bed"
done
```
Now separate bed files into upstream and downstream
```sh
mkdir -p tmp

FLANKS="/private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom/work/chm13v2.0.labels.as_hor.all_chrs.50kb_flanks.bed"
BED_DIR="/private/groups/patenlab/mira/centrolign/analysis/mutation_rate_comparison_to_flanks/per_sample_maf_beds"

# Split CHM13 flank BED into upstream (1st entry per chrom) and downstream (2nd entry per chrom)
awk '!seen[$1]++' "$FLANKS" > ./tmp/upstream_flanks.bed
awk 'seen[$1]++' "$FLANKS" > ./tmp/downstream_flanks.bed

# Split each sample mrg.bed into upstream and downstream
for bed in "$BED_DIR"/*.mrg.bed; do
    sample=$(basename "$bed" .mrg.bed)
    bedtools intersect -a "$bed" -b ./tmp/upstream_flanks.bed   > "${BED_DIR}/${sample}.upstream.bed"
    bedtools intersect -a "$bed" -b ./tmp/downstream_flanks.bed > "${BED_DIR}/${sample}.downstream.bed"
done
```

Get all pairs intersections of the up and downstream aligned beds per chrom
```sh
#!/bin/bash
#SBATCH --job-name=bedtools_intersect
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/subset_vcf_%x.%j.log
#SBATCH --time=12:00:00
#SBATCH --array=[0-23]%24

BED_DIR="/private/groups/patenlab/mira/centrolign/analysis/mutation_rate_comparison_to_flanks/per_sample_maf_beds"
OUT_DIR="/private/groups/patenlab/mira/centrolign/analysis/mutation_rate_comparison_to_flanks/pairwise_shared_aligned_bases"

source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate base

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")
chr=${chromosomes[$SLURM_ARRAY_TASK_ID]}

mkdir -p "$OUT_DIR"
OUT="$OUT_DIR/${chr}_pairwise_shared_aligned_bases.tsv"

echo -e "hap1\thap2\tflank\tshared_bases" > "$OUT"

haps=($(ls "$BED_DIR"/*.upstream.bed | xargs -I{} basename {} .upstream.bed))

n=${#haps[@]}
for ((i=0; i<n; i++)); do
    for ((j=i+1; j<n; j++)); do
        hap1="${haps[$i]}"
        hap2="${haps[$j]}"
        for flank in upstream downstream; do
            shared=$(bedtools intersect \
                -a "$BED_DIR/${hap1}.${flank}.bed" \
                -b "$BED_DIR/${hap2}.${flank}.bed" \
                | awk -v chr="$chr" '$1 == chr {sum += $3 - $2} END {print sum+0}')
            echo -e "${hap1}\t${hap2}\t${flank}\t${shared}"
        done
    done
done >> "$OUT"

```

#### Subsetting the MC vcf to 50kb up and downstream of the CHM13 arrays

```sh
#!/bin/bash
#SBATCH --job-name=subset_MC_vcf
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/subset_vcf_%x.%j.log
#SBATCH --time=1:00:00
#SBATCH --array=[0-23]%24

VCF="/private/groups/patenlab/mira/centrolign/analysis/mutation_rate_comparison_to_flanks/hprc-v2.0-mc-chm13.wave.vcf.gz"
BED_DIR="/private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom/work"
SAMPLE_DIR="/private/groups/migalab/juklucas/censat_regions/active_arrays"
OUTDIR="/private/groups/patenlab/mira/centrolign/analysis/mutation_rate_comparison_to_flanks/flanking_50kb"
FLANK=50000

mkdir -p ${OUTDIR}

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

chr=${chromosomes[$SLURM_ARRAY_TASK_ID]}

BED="${BED_DIR}/chm13v2.0.labels.as_hor.${chr}.active.mrg.bed"
SAMPLE_TSV="${SAMPLE_DIR}/asat_arrays_${chr}.tsv"

if [ ! -f "${BED}" ]; then
    echo "Skipping ${chr}: bed file not found"
    continue
fi
if [ ! -f "${SAMPLE_TSV}" ]; then
    echo "Skipping ${chr}: sample TSV not found"
    continue
fi

# Get alpha satellite boundaries
ASAT_START=$(awk 'NR>0 {print $2}' ${BED} | sort -n | head -1)
ASAT_END=$(awk 'NR>0 {print $3}' ${BED} | sort -n | tail -1)

LEFT_START=$(( ASAT_START - FLANK > 0 ? ASAT_START - FLANK : 0 ))
LEFT_END=${ASAT_START}
RIGHT_START=${ASAT_END}
RIGHT_END=$(( ASAT_END + FLANK ))

# Extract sample list from TSV (skip header, get unique sample_ids)
SAMPLE_FILE="${OUTDIR}/${chr}_samples.txt"
awk -F'\t' 'NR>1 {print $1}' ${SAMPLE_TSV} | sort -u > ${SAMPLE_FILE}

N_SAMPLES=$(wc -l < ${SAMPLE_FILE})
echo "${chr}: ${N_SAMPLES} samples, upstream ${LEFT_START}-${LEFT_END}, downstream ${RIGHT_START}-${RIGHT_END}"

bcftools view \
    -r "${chr}:${LEFT_START}-${LEFT_END},${chr}:${RIGHT_START}-${RIGHT_END}" \
    -S ${SAMPLE_FILE} \
    --force-samples \
    -O z -o "${OUTDIR}/${chr}_flanking_50kb.vcf.gz" \
    ${VCF}

bcftools index "${OUTDIR}/${chr}_flanking_50kb.vcf.gz"
```

```sh
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY; do
    VCF="/private/groups/patenlab/mira/centrolign/analysis/mutation_rate_comparison_to_flanks/flanking_50kb/${chr}_flanking_50kb.vcf.gz"
    [ ! -f "${VCF}" ] && continue

    TOTAL=$(bcftools query -l ${VCF} | wc -l)

    # Samples with any missing genotype
    N_MISSING=$(bcftools query -f '[%SAMPLE\t%GT\n]' ${VCF} | \
        awk '$2 ~ /\./' | cut -f1 | sort -u | wc -l)

    N_COMPLETE=$(( TOTAL - N_MISSING ))
    echo "${chr}: ${N_COMPLETE} / ${TOTAL} samples with no missing genotypes"
done
```

Remaining analysis done in notebooks/mutation_rate_comparison_to_flanks.ipynb
