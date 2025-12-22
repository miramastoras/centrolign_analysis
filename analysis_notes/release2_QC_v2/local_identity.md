### Intersecting centrolign pairwise variant calls with local identity annotations

Get list of SV beds
```sh
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

# master list of chr,SV_bed,
for chr in "${chromosomes[@]}";
  do grep "release 2 QC v2" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/centrolign_results.csv | grep -v "chr20," | grep -E "${chr},|${chr}_" | cut -f1,10 -d","  | while IFS=',' read -r subgroup cigar ; do
    find /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/${chr}/SV_beds/${subgroup}/ -type f | while read -r filepath; do
      echo "${chr},${subgroup},${filepath}"
    done
  done
done >> /private/groups/patenlab/mira/centrolign/analysis/local_identity/SVs/all_SV_beds.txt
```

```sh
#!/bin/bash
#SBATCH --job-name=local_identity_intersect
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/call_SVs_%x.%j.log
#SBATCH --time=1:00:00
#SBATCH --array=[1]%1

# activate environment for bedtools
source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate base

SV_BED_PATHS=$1 # csv file where each line is a path to the SV bedPE file for each pair

CHR=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$SV_BED_PATHS" | cut -f1 -d",")
SUBGROUP=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$SV_BED_PATHS" | cut -f1 -d",")
SV_BED=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$SV_BED_PATHS" | cut -f3 -d",")
REF_SMP=`echo $SV_BED | basename $SV_BED | cut -f1 -d"_"`
QUERY_SMP=`echo $SV_BED | basename $SV_BED | cut -f2 -d"_" | cut -f1-2 -d"."`

# index file for local identity beds
LOCAL_ID_IDX=/private/groups/migalab/juklucas/censat_paper/local_identity/2025_11_13/${CHR}/${CHR}_local_id_launch_outputs.csv

OUTDIR=/private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/local_identity/${SUBGROUP}/
mkdir -p $OUTDIR

echo $SV_BED
echo $REF_SMP
echo $QUERY_SMP

# convert bedPE into bed file
LOCAL_FOLDER=/data/tmp/$(whoami)/${REF_SMP}_${QUERY_SMP}_tmp/
mkdir -p ${LOCAL_FOLDER}

cut -f1-3,7,8 $SV_BED > ${LOCAL_FOLDER}/all_SVs.bed
cut -f4-8 $SV_BED >> ${LOCAL_FOLDER}/all_SVs.bed

# sep by triangle, trap, parallelogram
awk '$5 <= 0.1' ${LOCAL_FOLDER}/all_SVs.bed |  awk '$5 > 0' > ${LOCAL_FOLDER}/parallelograms.bed
awk '$5 > 0.1' ${LOCAL_FOLDER}/all_SVs.bed > ${LOCAL_FOLDER}/trapezoids.bed
awk '$5 == -1' ${LOCAL_FOLDER}/all_SVs.bed > ${LOCAL_FOLDER}/triangles.bed

# Get local identity bed file paths for each sample
REF_ID=`echo $REF_SMP | sed 's/\./#/g'`
QUERY_ID=`echo $QUERY_SMP | sed 's/\./#/g'`

REF_LOCAL_ID=`grep $REF_ID $LOCAL_ID_IDX | cut -f11 -d","`
QUERY_LOCAL_ID=`grep $REF_ID $LOCAL_ID_IDX | cut -f11 -d","`

# concatenate local identity tracks for both samples, convert contig name to match SV bed file
awk -v v="$REF_SMP" 'BEGIN { OFS="\t" } {$1=v; print}' ${REF_LOCAL_ID} > ${LOCAL_FOLDER}/local_identity_combined.bed
awk -v v="$QUERY_SMP" 'BEGIN { OFS="\t" } {$1=v; print}' ${QUERY_LOCAL_ID} >> ${LOCAL_FOLDER}/local_identity_combined.bed

# bedtools coverage to produce number of bases of each variant type overlapping each local identity window
bedtools coverage \
  -a ${LOCAL_FOLDER}/local_identity_combined.bed \
  -b ${LOCAL_FOLDER}/parallelograms.bed \
  > ${LOCAL_FOLDER}/local_iden_par.bed

bedtools coverage \
  -a ${LOCAL_FOLDER}/local_identity_combined.bed \
  -b ${LOCAL_FOLDER}/trapezoids.bed \
  > ${LOCAL_FOLDER}/local_iden_traps.bed

bedtools coverage \
  -a ${LOCAL_FOLDER}/local_identity_combined.bed \
  -b ${LOCAL_FOLDER}/triangles.bed \
  > ${LOCAL_FOLDER}/local_iden_tri.bed

cut -f1-9,11 ${LOCAL_FOLDER}/local_iden_par.bed > ${LOCAL_FOLDER}/par_counts.bed
cut -f11 ${LOCAL_FOLDER}/local_iden_traps.bed > ${LOCAL_FOLDER}/traps_counts.bed
cut -f11 ${LOCAL_FOLDER}/local_iden_tri.bed > ${LOCAL_FOLDER}/tri_counts.bed

paste -d"\t" ${LOCAL_FOLDER}/par_counts.bed ${LOCAL_FOLDER}/traps_counts.bed ${LOCAL_FOLDER}/tri_counts.bed > ${OUTDIR}/${REF_SMP}_${QUERY_SMP}.local_identity_SV_counts.bed

# clean up
rm -rf $LOCAL_FOLDER
```

```sh
cd /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/local_identity/
mkdir -p logs

# submit script to get count of bases within each local identity window for each sample pair
# output bed files will have # of bases from parallelograms, trapezoids, triangles as last 3 rows in that order
split -l 50000 -d -a 3 /private/groups/patenlab/mira/centrolign/analysis/local_identity/SVs/all_SV_beds.txt input_bed_lists/all_SV_beds_part_

sbatch --array=[1-50000]%150 \
    local_identity_intersect.sh \
    /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/local_identity/input_bed_lists/all_SV_beds_part_000

sbatch --array=[1-50000]%150 \
    local_identity_intersect.sh \
    /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/local_identity/input_bed_lists/all_SV_beds_part_001

sbatch --array=[1-50000]%150 \
    local_identity_intersect.sh \
    /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/local_identity/input_bed_lists/all_SV_beds_part_002

sbatch --array=[1-50000]%150 \
    local_identity_intersect.sh \
    /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/local_identity/input_bed_lists/all_SV_beds_part_003

sbatch --array=[1-50000]%150 \
    local_identity_intersect.sh \
    /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/local_identity/input_bed_lists/all_SV_beds_part_004
#!/bin/sh
sbatch --array=[1-50000]%150 \
    local_identity_intersect.sh \
    /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/local_identity/input_bed_lists/all_SV_beds_part_005

#!/bin/sh
sbatch --array=[1-50000]%150 \
    local_identity_intersect.sh \
    /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/local_identity/input_bed_lists/all_SV_beds_part_006

#!/bin/sh
sbatch --array=[1-50000]%150 \
    local_identity_intersect.sh \
    /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/local_identity/input_bed_lists/all_SV_beds_part_007

sbatch --array=[1-50000]%150 \
    local_identity_intersect.sh \
    /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/local_identity/input_bed_lists/all_SV_beds_part_008

#!/bin/sh
sbatch --array=[1-50000]%150 \
    local_identity_intersect.sh \
    /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/local_identity/input_bed_lists/all_SV_beds_part_009

sbatch --array=[1-36539]%150 \
    local_identity_intersect.sh \
    /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/local_identity/input_bed_lists/all_SV_beds_part_010
```
