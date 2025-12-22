### Checking alignability with respect to distance from CDR

#### Dist < 0.4

For all chromosomes, subgroups, call alignability windows for sample pairs < 0.4 distance
```sh
#!/bin/bash
#SBATCH --job-name=CDR_alignability_dist
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --array=[0]%24
#SBATCH --time=12:00:00

source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate ipynb

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

chr=${chromosomes[$SLURM_ARRAY_TASK_ID]}

mkdir -p /private/groups/patenlab/mira/centrolign/analysis/alignability/windowed_align_beds_50bp_d0.4/${chr}/

# call alignability windows
grep "release 2 QC v2" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/centrolign_results.csv | grep -v "chr20," | grep -E "${chr},|${chr}_" | cut -f1,8,10 -d","  | while IFS=',' read -r subgroup fasta cigar ; do
  time python /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/notebooks/call_alignability_windows.py \
  --cigar-path ${cigar} \
  --fai ${fasta}.fai \
  --distance-csv /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/${chr}_r2_QC_v2_centrolign_pairwise_distance.csv \
  --distance-threshold 0.4 \
  --window-size 50 \
  --outdir /private/groups/patenlab/mira/centrolign/analysis/alignability/windowed_align_beds_50bp_d0.4/${chr}/
done

# Convert chr11 array coordinates to assembly coordinates
mkdir -p /private/groups/patenlab/mira/centrolign/analysis/alignability/windowed_align_beds_50bp_d0.4_asm_coords/${chr}/

time python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/convert_SV_bed_to_asm_coords.py \
  --format bed \
  -s /private/groups/patenlab/mira/centrolign/analysis/alignability/windowed_align_beds_50bp_d0.4/${chr}/ \
  -a /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/ \
  -c ${chr} \
  -o /private/groups/patenlab/mira/centrolign/analysis/alignability/windowed_align_beds_50bp_d0.4_asm_coords/${chr}/

# Use bedtools closest to calculate distance from nearest CDR
conda activate base

cd /private/groups/patenlab/mira/centrolign/analysis/alignability/windowed_align_beds_50bp_d0.4_asm_coords/${chr}/

mkdir -p /private/groups/patenlab/mira/centrolign/analysis/alignability/windowed_align_beds_50bp_d0.4_asm_coords/${chr}/CDR_distance/

ls *.bed | while read -r line; do
  bedName=$(basename "$line" .bed)

  bedtools closest \
    -a "$line" \
    -b /private/groups/migalab/jmmenend/HPRC/cenSatProject/CDR_data/${chr}/${chr}.centrodip.bed \
    -D a -t first \
    > /private/groups/patenlab/mira/centrolign/analysis/alignability/windowed_align_beds_50bp_d0.4_asm_coords/${chr}/CDR_distance/${bedName}.CDR_dist.bed
done
```

#### Dist < 0.2

For all chromosomes, subgroups, call alignability windows for sample pairs < 0.2 distance
```sh
#!/bin/bash
#SBATCH --job-name=CDR_alignability_dist
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --array=[0]%24
#SBATCH --time=12:00:00

source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate ipynb

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

chr=${chromosomes[$SLURM_ARRAY_TASK_ID]}

mkdir -p /private/groups/patenlab/mira/centrolign/analysis/alignability/windowed_align_beds_50bp_d0.2/${chr}/

# call alignability windows
grep "release 2 QC v2" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/centrolign_results.csv | grep -v "chr20," | grep -E "${chr},|${chr}_" | cut -f1,8,10 -d","  | while IFS=',' read -r subgroup fasta cigar ; do
  time python /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/notebooks/call_alignability_windows.py \
  --cigar-path ${cigar} \
  --fai ${fasta}.fai \
  --distance-csv /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/${chr}_r2_QC_v2_centrolign_pairwise_distance.csv \
  --distance-threshold 0.2 \
  --window-size 50 \
  --outdir /private/groups/patenlab/mira/centrolign/analysis/alignability/windowed_align_beds_50bp_d0.2/${chr}/
done

# Convert chr11 array coordinates to assembly coordinates
mkdir -p /private/groups/patenlab/mira/centrolign/analysis/alignability/windowed_align_beds_50bp_d0.2_asm_coords/${chr}/

time python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/convert_SV_bed_to_asm_coords.py \
  --format bed \
  -s /private/groups/patenlab/mira/centrolign/analysis/alignability/windowed_align_beds_50bp_d0.2/${chr}/ \
  -a /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/ \
  -c ${chr} \
  -o /private/groups/patenlab/mira/centrolign/analysis/alignability/windowed_align_beds_50bp_d0.2_asm_coords/${chr}/

# Use bedtools closest to calculate distance from nearest CDR
conda activate base

cd /private/groups/patenlab/mira/centrolign/analysis/alignability/windowed_align_beds_50bp_d0.2_asm_coords/${chr}/

mkdir -p /private/groups/patenlab/mira/centrolign/analysis/alignability/windowed_align_beds_50bp_d0.2_asm_coords/${chr}/CDR_distance/

ls *.bed | while read -r line; do
  bedName=$(basename "$line" .bed)

  bedtools closest \
    -a "$line" \
    -b /private/groups/migalab/jmmenend/HPRC/cenSatProject/CDR_data/${chr}/${chr}.centrodip.bed \
    -D a -t first \
    > /private/groups/patenlab/mira/centrolign/analysis/alignability/windowed_align_beds_50bp_d0.2_asm_coords/${chr}/CDR_distance/${bedName}.CDR_dist.bed
done
```
