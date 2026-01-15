### Linear regression of CDR distance, local homogeneity, chromosome

Get distance from CDR for local identity windows
```sh
#!/bin/bash
#SBATCH --job-name=local_identity_CDR_dist
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x.%j.log
#SBATCH --time=1:00:00
#SBATCH --array=[0-23]%24

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

CHR=${chromosomes[$SLURM_ARRAY_TASK_ID]}

# activate environment for bedtools
source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate base

# index file for local identity beds
LOCAL_ID_IDX=/private/groups/migalab/juklucas/censat_paper/local_identity/2025_11_13/${CHR}/${CHR}_local_id_launch_outputs.csv

LOCAL_FOLDER=/data/tmp/$(whoami)/${CHR}_cdr_dist_local_id/
OUTDIR=/private/groups/patenlab/mira/centrolign/analysis/CDR_variant_regression/local_id_CDR_dist/${CHR}
mkdir -p ${LOCAL_FOLDER}
mkdir -p ${OUTDIR}

# For each sample, get distance to CDR for local identity bed
grep -v "chm13" $LOCAL_ID_IDX | grep -v "hg002v1.1" | cut -f4,11 -d"," | tail -n +2 | while IFS=',' read -r contig local_id_bed ; do

  SAMPLE=`echo $contig | cut -f1-2 -d"#" | sed 's/#/./g'`

  # convert local identity bed to assembly coordinates
  ASAT_START=`grep "${CHR}$" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/${SAMPLE}_asat_arrays.bed | cut -f2`

  awk -F'\t' -v OFS='\t' -v n=$ASAT_START '{ $2 += n; $3 += n }1' ${local_id_bed} > ${LOCAL_FOLDER}/${SAMPLE}_local_id_asat_coords.bed

  # get distance to CDR
  bedtools closest \
    -a ${LOCAL_FOLDER}/${SAMPLE}_local_id_asat_coords.bed \
    -b /private/groups/migalab/jmmenend/HPRC/cenSatProject/CDR_data/${CHR}/${CHR}.centrodip.bed \
    -D a -t first \
    > ${LOCAL_FOLDER}/${SAMPLE}.local_id_CDR_dist.bed

  # just save relevant columns
  #awk -F'\t' 'BEGIN {OFS="\t"} {print $1, $2, $3, $5, $19}' ${LOCAL_FOLDER}/${SAMPLE}.local_id_CDR_dist.bed | sed "s/$contig/$SAMPLE/g" > ${OUTDIR}/${SAMPLE}.local_id_CDR_dist.bed

  awk -F'\t' 'BEGIN {OFS="\t"} {print $1, $2, $3, $5, $19}' ${LOCAL_FOLDER}/${SAMPLE}.local_id_CDR_dist.bed > ${OUTDIR}/${SAMPLE}.local_id_CDR_dist.bed

done

rm -rf $LOCAL_FOLDER
```

Calculate SNV rate within local identity windows

```sh
#!/bin/bash
#SBATCH --job-name=CDR_dist_snv_windows_filt
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/call_SVs_%x.%j.log
#SBATCH --time=12:00:00
#SBATCH --array=[1]%24

# activate environment for bedtools
source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate base

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

CHR=${chromosomes[$SLURM_ARRAY_TASK_ID]}

LOCAL_FOLDER=/data/tmp/$(whoami)/${CHR}_SNVs_tmps/
mkdir -p ${LOCAL_FOLDER}
OUTDIR=/private/groups/patenlab/mira/centrolign/analysis/CDR_variant_regression/SNVs_pairwise_raw/${CHR}
mkdir -p ${OUTDIR}

cd /private/groups/patenlab/mira/centrolign/analysis/SNVs_pairwise_asm_coords/${CHR}/

ls *.bed | head -n 1 | while read -r bed ; do
        # for each bed get ref and query sample ID
        REF_SMP=`echo $bed | basename $bed | cut -f1 -d"_"`
        QUERY_SMP=`echo $bed | basename $bed | cut -f2 -d"_" | cut -f1-2 -d"."`

        # get asm contig name
        SMP1_contig=`head -n1 $bed | cut -f1`
        SMP2_contig=`tail -n1 $bed | cut -f1`

        # get file identifier - remove extension
        FILE_ID=${bed%.*}

        # concatenate window files for both samples
        cat /private/groups/patenlab/mira/centrolign/analysis/CDR_variant_regression/local_id_CDR_dist/${CHR}/${REF_SMP}.local_id_CDR_dist.bed /private/groups/patenlab/mira/centrolign/analysis/CDR_variant_regression/local_id_CDR_dist/${CHR}/${QUERY_SMP}.local_id_CDR_dist.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.local_id_CDR_dist.bed

        # Use bedtools coverage to calculate number of SNV bases overlapping windows
        bedtools coverage \
          -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.local_id_CDR_dist.bed \
          -b ${bed} \
          > ${OUTDIR}/${FILE_ID}.all_snvs.local_id_CDR_dist.bed
    done

rm -rf ${LOCAL_FOLDER}
```

Calculate short indel triangle rate within local identity windows

```sh
#!/bin/bash
#SBATCH --job-name=CDR_dist_short_indel_windows
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x.%j.log
#SBATCH --time=12:00:00
#SBATCH --array=[1]%24

# activate environment for bedtools
source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate base

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

CHR=${chromosomes[$SLURM_ARRAY_TASK_ID]}

LOCAL_FOLDER=/data/tmp/$(whoami)/${CHR}_short_indels_tmps/
mkdir -p ${LOCAL_FOLDER}
OUTDIR=/private/groups/patenlab/mira/centrolign/analysis/CDR_variant_regression/short_indels_triangles/${CHR}
mkdir -p ${OUTDIR}

grep "release 2 QC v2" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/centrolign_results.csv | grep -v "chr20," | grep -E "${CHR},|${CHR}_" | cut -f1,8,10 -d","  | head -n 1 | while IFS=',' read -r subgroup fasta cigar ; do

    LOCAL_FOLDER=/data/tmp/$(whoami)/${CHR}_${subgroup}_short_indels_tmp/
    mkdir -p ${LOCAL_FOLDER}
    mkdir -p /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/short_indels_pairwise/${WINDOWSIZE}/${CHR}/

    cd /private/groups/patenlab/mira/centrolign/analysis/short_indels_pairwise_asm_coords/${CHR}/short_indel_beds/${subgroup}/
    ls *.bed | head -n 1 | while read -r bed ; do
        # for each bed get ref and query sample ID
        REF_SMP=`echo $bed | basename $bed | cut -f1 -d"_"`
        QUERY_SMP=`echo $bed | basename $bed | cut -f2 -d"_" | cut -f1-2 -d"."`

        # convert bedPE into bed file, containing SV calls for both samples
        cut -f1-3,7,8 $bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_indels.bed
        cut -f4-8 $bed >> ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_indels.bed

        SMP1_contig=`cut -f1 $bed | sort | uniq`
        SMP2_contig=`cut -f4 $bed | sort | uniq`

        # sep SVs by triangle, trap, parallelogram
        awk '$5 <= 0.1' ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_indels.bed |  awk '$5 >= 0' > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.parallelograms.bed
        awk '$5 > 0.1' ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_indels.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.trapezoids.bed
        awk '$5 == -1' ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_indels.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.triangles.bed

        # concatenate window files for both samples,
        cat /private/groups/patenlab/mira/centrolign/analysis/CDR_variant_regression/local_id_CDR_dist/${CHR}/${REF_SMP}.local_id_CDR_dist.bed /private/groups/patenlab/mira/centrolign/analysis/CDR_variant_regression/local_id_CDR_dist/${CHR}/${QUERY_SMP}.local_id_CDR_dist.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.local_id_CDR_dist.bed

        # Use bedtools coverage to calculate number of SV bases overlapping windows
        bedtools coverage \
          -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.local_id_CDR_dist.bed \
          -b ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.parallelograms.bed \
          > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.local_id_CDR.par.bed

        bedtools coverage \
          -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.local_id_CDR_dist.bed \
          -b ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.trapezoids.bed \
          > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.local_id_CDR.trap.bed

        bedtools coverage \
          -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.local_id_CDR_dist.bed \
          -b ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.triangles.bed \
          > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.local_id_CDR.tri.bed

        # combine into single file
        cut -f1-6,9 ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.local_id_CDR.par.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.par_counts.bed
        cut -f6,9 ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.local_id_CDR.trap.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.trap_counts.bed
        cut -f6,9 ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.local_id_CDR.tri.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.tri_counts.bed

        paste -d"\t" ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.par_counts.bed ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.trap_counts.bed ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.tri_counts.bed > ${OUTDIR}/${REF_SMP}_${QUERY_SMP}.all_short_indel_counts.bed

    done

    rm -rf ${LOCAL_FOLDER}
done
```

Calculate SV triangle rate within local identity windows

```sh
#!/bin/bash
#SBATCH --job-name=CDR_dist_SV_windows
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x.%j.log
#SBATCH --time=12:00:00
#SBATCH --array=[0-23]%24

# activate environment for bedtools
source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate base

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

CHR=${chromosomes[$SLURM_ARRAY_TASK_ID]}

LOCAL_FOLDER=/data/tmp/$(whoami)/${CHR}_SVs_tmps/
mkdir -p ${LOCAL_FOLDER}
OUTDIR=/private/groups/patenlab/mira/centrolign/analysis/CDR_variant_regression/SVs_pairwise_all/${CHR}
mkdir -p ${OUTDIR}

grep "release 2 QC v2" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/centrolign_results.csv | grep -v "chr20," | grep -E "${CHR},|${CHR}_" | cut -f1,8,10 -d","  | head -n 1 | while IFS=',' read -r subgroup fasta cigar ; do

    LOCAL_FOLDER=/data/tmp/$(whoami)/${CHR}_${subgroup}_SVs_tmp/
    mkdir -p ${LOCAL_FOLDER}

    cd /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise_asm_coords/${CHR}/SV_beds/${subgroup}/
    ls *.bed | head -n 1 | while read -r bed ; do
        # for each bed get ref and query sample ID
        REF_SMP=`echo $bed | basename $bed | cut -f1 -d"_"`
        QUERY_SMP=`echo $bed | basename $bed | cut -f2 -d"_" | cut -f1-2 -d"."`

        # convert bedPE into bed file, containing SV calls for both samples
        cut -f1-3,7,8 $bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_SVs.bed
        cut -f4-8 $bed >> ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_SVs.bed

        SMP1_contig=`cut -f1 $bed | sort | uniq`
        SMP2_contig=`cut -f4 $bed | sort | uniq`

        # sep SVs by triangle, trap, parallelogram
        awk '$5 <= 0.1' ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_SVs.bed |  awk '$5 >= 0' > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.parallelograms.bed
        awk '$5 > 0.1' ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_SVs.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.trapezoids.bed
        awk '$5 == -1' ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_SVs.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.triangles.bed

        # concatenate window files for both samples,
        cat /private/groups/patenlab/mira/centrolign/analysis/CDR_variant_regression/local_id_CDR_dist/${CHR}/${REF_SMP}.local_id_CDR_dist.bed /private/groups/patenlab/mira/centrolign/analysis/CDR_variant_regression/local_id_CDR_dist/${CHR}/${QUERY_SMP}.local_id_CDR_dist.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.local_id_CDR_dist.bed

        # Use bedtools coverage to calculate number of SV bases overlapping windows
        bedtools coverage \
          -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.local_id_CDR_dist.bed \
          -b ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.parallelograms.bed \
          > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.local_id_CDR.par.bed

        bedtools coverage \
          -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.local_id_CDR_dist.bed \
          -b ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.trapezoids.bed \
          > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.local_id_CDR.trap.bed

        bedtools coverage \
          -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.local_id_CDR_dist.bed \
          -b ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.triangles.bed \
          > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.local_id_CDR.tri.bed

        # combine into single file
        cut -f1-6,9 ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.local_id_CDR.par.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.par_counts.bed
        cut -f6,9 ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.local_id_CDR.trap.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.trap_counts.bed
        cut -f6,9 ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.local_id_CDR.tri.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.tri_counts.bed

        paste -d"\t" ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.par_counts.bed ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.trap_counts.bed ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.tri_counts.bed > ${OUTDIR}/${REF_SMP}_${QUERY_SMP}.all_SV_counts.bed

    done

    rm -rf ${LOCAL_FOLDER}
done
```

### Use CDR midpoint

Convert CDR bed files to CDR midpoints
```sh
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

for chr in "${chromosomes[@]}";
  do
  awk 'BEGIN{OFS="\t"} {mid=int(($2+$3)/2); print $1, mid, mid+1, $4, $5, $6, $7, $8, $9}' /private/groups/migalab/jmmenend/HPRC/cenSatProject/CDR_data/${chr}/${chr}.centrodip.bed > /private/groups/patenlab/mira/centrolign/annotations/centrodip_midpoints_01132025/${chr}.centrodip.bed
done
```

Get distance from CDR for local identity windows
```sh
#!/bin/bash
#SBATCH --job-name=local_identity_CDR_dist
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x.%j.log
#SBATCH --time=1:00:00
#SBATCH --array=[0-23]%24

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

CHR=${chromosomes[$SLURM_ARRAY_TASK_ID]}

# activate environment for bedtools
source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate base

# index file for local identity beds
LOCAL_ID_IDX=/private/groups/migalab/juklucas/censat_paper/local_identity/2025_11_13/${CHR}/${CHR}_local_id_launch_outputs.csv

LOCAL_FOLDER=/data/tmp/$(whoami)/${CHR}_cdr_dist_local_id/
OUTDIR=/private/groups/patenlab/mira/centrolign/analysis/CDR_variant_regression/local_id_CDR_dist_midpoint/${CHR}
mkdir -p ${LOCAL_FOLDER}
mkdir -p ${OUTDIR}

# For each sample, get distance to CDR for local identity bed
grep -v "chm13" $LOCAL_ID_IDX | grep -v "hg002v1.1" | cut -f4,11 -d"," | tail -n +2 | while IFS=',' read -r contig local_id_bed ; do

  SAMPLE=`echo $contig | cut -f1-2 -d"#" | sed 's/#/./g'`

  # convert local identity bed to assembly coordinates
  ASAT_START=`grep "${CHR}$" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/${SAMPLE}_asat_arrays.bed | cut -f2`

  awk -F'\t' -v OFS='\t' -v n=$ASAT_START '{ $2 += n; $3 += n }1' ${local_id_bed} > ${LOCAL_FOLDER}/${SAMPLE}_local_id_asat_coords.bed

  # get distance to CDR
  bedtools closest \
    -a ${LOCAL_FOLDER}/${SAMPLE}_local_id_asat_coords.bed \
    -b /private/groups/patenlab/mira/centrolign/annotations/centrodip_midpoints_01132025/${CHR}.centrodip.bed \
    -D a -t first \
    > ${LOCAL_FOLDER}/${SAMPLE}.local_id_CDR_dist.bed

  # just save relevant columns
  #awk -F'\t' 'BEGIN {OFS="\t"} {print $1, $2, $3, $5, $19}' ${LOCAL_FOLDER}/${SAMPLE}.local_id_CDR_dist.bed | sed "s/$contig/$SAMPLE/g" > ${OUTDIR}/${SAMPLE}.local_id_CDR_dist.bed

  awk -F'\t' 'BEGIN {OFS="\t"} {print $1, $2, $3, $5, $19}' ${LOCAL_FOLDER}/${SAMPLE}.local_id_CDR_dist.bed > ${OUTDIR}/${SAMPLE}.local_id_CDR_dist.bed

done

rm -rf $LOCAL_FOLDER
```

Calculate SNV rate within local identity windows

```sh
#!/bin/bash
#SBATCH --job-name=CDR_dist_snv_windows_filt
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/call_SVs_%x.%j.log
#SBATCH --time=12:00:00
#SBATCH --array=[0-23]%24

# activate environment for bedtools
source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate base

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

CHR=${chromosomes[$SLURM_ARRAY_TASK_ID]}

LOCAL_FOLDER=/data/tmp/$(whoami)/${CHR}_SNVs_tmps/
mkdir -p ${LOCAL_FOLDER}
OUTDIR=/private/groups/patenlab/mira/centrolign/analysis/CDR_variant_regression/SNVs_pairwise_raw_midpoint/${CHR}
mkdir -p ${OUTDIR}

cd /private/groups/patenlab/mira/centrolign/analysis/SNVs_pairwise_asm_coords/${CHR}/

ls *.bed | while read -r bed ; do
        # for each bed get ref and query sample ID
        REF_SMP=`echo $bed | basename $bed | cut -f1 -d"_"`
        QUERY_SMP=`echo $bed | basename $bed | cut -f2 -d"_" | cut -f1-2 -d"."`

        # get asm contig name
        SMP1_contig=`head -n1 $bed | cut -f1`
        SMP2_contig=`tail -n1 $bed | cut -f1`

        # get file identifier - remove extension
        FILE_ID=${bed%.*}

        # concatenate window files for both samples
        cat /private/groups/patenlab/mira/centrolign/analysis/CDR_variant_regression/local_id_CDR_dist_midpoint/${CHR}/${REF_SMP}.local_id_CDR_dist.bed /private/groups/patenlab/mira/centrolign/analysis/CDR_variant_regression/local_id_CDR_dist_midpoint/${CHR}/${QUERY_SMP}.local_id_CDR_dist.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.local_id_CDR_dist.bed

        # Use bedtools coverage to calculate number of SNV bases overlapping windows
        bedtools coverage \
          -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.local_id_CDR_dist.bed \
          -b ${bed} \
          > ${OUTDIR}/${FILE_ID}.all_snvs.local_id_CDR_dist.bed
    done

rm -rf ${LOCAL_FOLDER}
```
