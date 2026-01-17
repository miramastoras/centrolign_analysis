# This document contains analysis for plotting mutation rate as a function of distance to the CDR

## For all SVs, calculate their distance to CDR

Convert all SV beds into assembly coordinates
```sh
#!/bin/bash
#SBATCH --job-name=convert_SV_bed_to_asm_coords
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --array=[0-23]%24
#SBATCH --time=12:00:00

source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate ipynb

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

chr=${chromosomes[$SLURM_ARRAY_TASK_ID]}

grep "release 2 QC v2" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/centrolign_results.csv | grep -v "chr20," | grep -E "${chr},|${chr}_" | cut -f1,8,10 -d","  | while IFS=',' read -r subgroup fasta cigar ; do
  mkdir -p /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise_asm_coords/${chr}/SV_beds/${subgroup}/

  time python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/convert_SV_bed_to_asm_coords.py \
  --format bedpe \
  -s /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/${chr}/SV_beds/${subgroup}/ \
  -a /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/ \
  -c ${chr} \
  -o /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise_asm_coords/${chr}/SV_beds/${subgroup}/

done
```

Take SV bed file and calculate distance to CDR
```sh
#!/bin/bash
#SBATCH --job-name=CDR_dist_SVs
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/call_SVs_%x.%j.log
#SBATCH --time=12:00:00
#SBATCH --array=23

# activate environment for bedtools
source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate base

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

chr=${chromosomes[$SLURM_ARRAY_TASK_ID]}

# Use bedtools closest to calculate distance from the CDR
conda activate base

grep "release 2 QC v2" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/centrolign_results.csv | grep -v "chr20," | grep -E "${chr},|${chr}_" | cut -f1,8,10 -d","  | while IFS=',' read -r subgroup fasta cigar ; do

    LOCAL_FOLDER=/data/tmp/$(whoami)/${chr}_${subgroup}_tmp/
    mkdir -p ${LOCAL_FOLDER}

    cd /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise_asm_coords/${chr}/SV_beds/${subgroup}/

    mkdir -p /private/groups/patenlab/mira/centrolign/analysis/variant_dist_to_CDR/${chr}

    ls *.bed | while read -r bed ; do
        REF_SMP=`echo $bed | basename $bed | cut -f1 -d"_"`
        QUERY_SMP=`echo $bed | basename $bed | cut -f2 -d"_" | cut -f1-2 -d"."`

        # convert bedPE into bed file
        cut -f1-3,7,8 $bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_SVs.bed
        cut -f4-8 $bed >> ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_SVs.bed

        bedtools sort -i ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_SVs.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_SVs.srt.bed

        # Get distance to CDR for every SV call
        bedtools closest \
          -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_SVs.srt.bed \
          -b /private/groups/migalab/jmmenend/HPRC/cenSatProject/CDR_data/${chr}/${chr}.centrodip.bed \
          -D a -t first \
          > /private/groups/patenlab/mira/centrolign/analysis/variant_dist_to_CDR/${chr}/${REF_SMP}_${QUERY_SMP}.CDR_dist.bed
    done

    # clean up tmp beds
    rm -rf ${LOCAL_FOLDER}
done
```

## Starting over: using windows instead to calculate rates

For each array, make windows of 100 bp, 1kb, 10kb, 100kb.
Count number of SV bases overlapping those windows per sample pair
Get rate per window (# of SV bases / window size)
Use bedtools closest to assign each window to a distance from the CDR.
As with the alignability, bin the windows within distance from CDR and plot average rate per window.

Make windows across the asat arrays in asm coords for all samples

```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds
mkdir -p /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/per_smp_asat_beds_windows/10kb/
mkdir -p /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/per_smp_asat_beds_windows/100kb/

ls *.bed | while read -r bed ; do
  bedtools makewindows -b ${bed} -w 10000 > /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/per_smp_asat_beds_windows/10kb/${bed}
  bedtools makewindows -b ${bed} -w 100000 > /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/per_smp_asat_beds_windows/100kb/${bed}
done
```

Script to get count of SVs within each window, and then windows distance from CDR

```sh
#!/bin/bash
#SBATCH --job-name=CDR_dist_SVs_windows
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/call_SVs_%x.%j.log
#SBATCH --time=12:00:00
#SBATCH --array=[23]%24

# activate environment for bedtools
source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate base

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

chr=${chromosomes[$SLURM_ARRAY_TASK_ID]}

WINDOWSIZE="100kb"

grep "release 2 QC v2" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/centrolign_results.csv | grep -v "chr20," | grep -E "${chr},|${chr}_" | cut -f1,8,10 -d","  | while IFS=',' read -r subgroup fasta cigar ; do

    LOCAL_FOLDER=/data/tmp/$(whoami)/${chr}_${subgroup}_tmp/
    mkdir -p ${LOCAL_FOLDER}
    mkdir -p /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/SVs_induced_pairwise/${WINDOWSIZE}/${chr}/

    cd /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise_asm_coords/${chr}/SV_beds/${subgroup}/
    ls *.bed | while read -r bed ; do
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

        # concatenate window files for both samples, restricting just to contigs for this chrom
        grep $SMP1_contig /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/per_smp_asat_beds_windows/${WINDOWSIZE}/${REF_SMP}_asat_arrays.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.bed

        grep $SMP2_contig /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/per_smp_asat_beds_windows/${WINDOWSIZE}/${QUERY_SMP}_asat_arrays.bed >> ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.bed

        # Use bedtools coverage to calculate number of SV bases overlapping windows
        bedtools coverage \
          -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.bed \
          -b ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.parallelograms.bed \
          > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.par.bed

        bedtools coverage \
          -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.bed \
          -b ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.trapezoids.bed \
          > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.trap.bed

        bedtools coverage \
          -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.bed \
          -b ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.triangles.bed \
          > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.tri.bed

        cut -f1-3,5,7 ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.par.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.par_counts.bed
        cut -f5,7 ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.trap.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.trap_counts.bed
        cut -f5,7 ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.tri.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.tri_counts.bed

        paste -d"\t" ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.par_counts.bed ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.trap_counts.bed ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.tri_counts.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.all_SV_counts.bed

        # use bedtools closest to get distance from CDR
        bedtools sort -i ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.all_SV_counts.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.all_SV_counts.srt.bed

        # Get distance to CDR for every SV call
        bedtools closest \
          -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.all_SV_counts.srt.bed \
          -b /private/groups/migalab/jmmenend/HPRC/cenSatProject/CDR_data/${chr}/${chr}.centrodip.bed \
          -D a -t first \
          > /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/SVs_induced_pairwise/${WINDOWSIZE}/${chr}/${REF_SMP}_${QUERY_SMP}.CDR_dist.${WINDOWSIZE}.bed
    done

    rm -rf ${LOCAL_FOLDER}
done
```

### Make sure to get count of SVs / samples excluded because there is no CDR intersection - may need to exclude the entire sample pair?

### Run for short indels
```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds
mkdir -p /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/per_smp_asat_beds_windows/100bp/
mkdir -p /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/per_smp_asat_beds_windows/1kb/

ls *.bed | while read -r bed ; do
  bedtools makewindows -b ${bed} -w 100 > /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/per_smp_asat_beds_windows/100bp/${bed}
  bedtools makewindows -b ${bed} -w 1000 > /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/per_smp_asat_beds_windows/1kb/${bed}
done
```

Convert all SV beds into assembly coordinates
```sh
#!/bin/bash
#SBATCH --job-name=convert_short_indel_bed_to_asm_coords
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --array=[0-23]%24
#SBATCH --time=12:00:00

source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate ipynb

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

chr=${chromosomes[$SLURM_ARRAY_TASK_ID]}

grep "release 2 QC v2" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/centrolign_results.csv | grep -v "chr20," | grep -E "${chr},|${chr}_" | cut -f1,8,10 -d","  | while IFS=',' read -r subgroup fasta cigar ; do
  mkdir -p /private/groups/patenlab/mira/centrolign/analysis/short_indels_pairwise_asm_coords/${chr}/short_indel_beds/${subgroup}/

  time python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/convert_SV_bed_to_asm_coords.py \
  --format bedpe \
  -s /private/groups/patenlab/mira/centrolign/analysis/short_indels_pairwise/${chr}/short_indel_beds/${subgroup}/ \
  -a /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/ \
  -c ${chr} \
  -o /private/groups/patenlab/mira/centrolign/analysis/short_indels_pairwise_asm_coords/${chr}/short_indel_beds/${subgroup}/

done
```

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
#SBATCH --output=logs/call_SVs_%x.%j.log
#SBATCH --time=12:00:00
#SBATCH --array=[0-23]%24

# activate environment for bedtools
source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate base

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

chr=${chromosomes[$SLURM_ARRAY_TASK_ID]}

WINDOWSIZE="1kb"

grep "release 2 QC v2" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/centrolign_results.csv | grep -v "chr20," | grep -E "${chr},|${chr}_" | cut -f1,8,10 -d","  | while IFS=',' read -r subgroup fasta cigar ; do

    LOCAL_FOLDER=/data/tmp/$(whoami)/${chr}_${subgroup}_short_indels_tmp/
    mkdir -p ${LOCAL_FOLDER}
    mkdir -p /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/short_indels_pairwise/${WINDOWSIZE}/${chr}/

    cd /private/groups/patenlab/mira/centrolign/analysis/short_indels_pairwise_asm_coords/${chr}/short_indel_beds/${subgroup}/
    ls *.bed | while read -r bed ; do
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

        # concatenate window files for both samples, restricting just to contigs for this chrom
        grep $SMP1_contig /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/per_smp_asat_beds_windows/${WINDOWSIZE}/${REF_SMP}_asat_arrays.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.bed

        grep $SMP2_contig /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/per_smp_asat_beds_windows/${WINDOWSIZE}/${QUERY_SMP}_asat_arrays.bed >> ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.bed

        # Use bedtools coverage to calculate number of SV bases overlapping windows
        bedtools coverage \
          -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.bed \
          -b ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.parallelograms.bed \
          > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.par.bed

        bedtools coverage \
          -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.bed \
          -b ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.trapezoids.bed \
          > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.trap.bed

        bedtools coverage \
          -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.bed \
          -b ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.triangles.bed \
          > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.tri.bed

        cut -f1-3,5,7 ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.par.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.par_counts.bed
        cut -f5,7 ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.trap.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.trap_counts.bed
        cut -f5,7 ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.tri.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.tri_counts.bed

        paste -d"\t" ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.par_counts.bed ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.trap_counts.bed ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.tri_counts.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.all_short_indel_counts.bed

        # use bedtools closest to get distance from CDR
        bedtools sort -i ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.all_short_indel_counts.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.all_short_indel_counts.srt.bed

        # Get distance to CDR for every SV call
        bedtools closest \
          -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.all_short_indel_counts.srt.bed \
          -b /private/groups/migalab/jmmenend/HPRC/cenSatProject/CDR_data/${chr}/${chr}.centrodip.bed \
          -D a -t first \
          > /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/short_indels_pairwise/${WINDOWSIZE}/${chr}/${REF_SMP}_${QUERY_SMP}.CDR_dist.${WINDOWSIZE}.bed
    done

    rm -rf ${LOCAL_FOLDER}
done
```

Run for SNVs

### Run for SNVs - raw

Convert Julian's SNV csv files into bed files in asm coordinates
- reran this correctig my script convert_SNV_csv_to_bed.py because Julian's files weren't bed coordinates, everything was shifted by 1
```sh
#!/bin/bash
#SBATCH --job-name=convert_SNV_bed_to_asm_coords
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --array=[0-23]%24
#SBATCH --time=12:00:00

source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate ipynb

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

chr=${chromosomes[$SLURM_ARRAY_TASK_ID]}

mkdir -p /private/groups/patenlab/mira/centrolign/analysis/SNVs_pairwise_asm_coords/${chr}/

time python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/convert_SNV_csv_to_bed.py \
  -i /private/groups/migalab/juklucas/centrolign/variant_calling/rates/${chr}/snv_calls/raw/ \
  -a /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/ \
  -c ${chr} \
  -o /private/groups/patenlab/mira/centrolign/analysis/SNVs_pairwise_asm_coords/${chr}/
```

```sh
#!/bin/bash
#SBATCH --job-name=CDR_dist_snv_windows
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

chr=${chromosomes[$SLURM_ARRAY_TASK_ID]}

WINDOWSIZE="1kb"

LOCAL_FOLDER=/data/tmp/$(whoami)/${chr}_SNVs_tmps/
mkdir -p ${LOCAL_FOLDER}
mkdir -p /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/SNVs_pairwise_raw/${WINDOWSIZE}/${chr}

cd /private/groups/patenlab/mira/centrolign/analysis/SNVs_pairwise_asm_coords/${chr}/

ls *.bed | while read -r bed ; do
        # for each bed get ref and query sample ID
        REF_SMP=`echo $bed | basename $bed | cut -f1 -d"_"`
        QUERY_SMP=`echo $bed | basename $bed | cut -f2 -d"_" | cut -f1-2 -d"."`

        # get asm contig name
        SMP1_contig=`head -n1 $bed | cut -f1`
        SMP2_contig=`tail -n1 $bed | cut -f1`

        # concatenate window files for both samples, restricting just to contigs for this chrom
        grep $SMP1_contig /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/per_smp_asat_beds_windows/${WINDOWSIZE}/${REF_SMP}_asat_arrays.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.bed

        grep $SMP2_contig /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/per_smp_asat_beds_windows/${WINDOWSIZE}/${QUERY_SMP}_asat_arrays.bed >> ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.bed

        # Use bedtools coverage to calculate number of SNV bases overlapping windows
        bedtools coverage \
          -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.bed \
          -b ${bed} \
          > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_snvs.w${WINDOWSIZE}.bed

        # use bedtools closest to get distance from CDR
        bedtools sort -i ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_snvs.w${WINDOWSIZE}.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_snvs.w${WINDOWSIZE}.srt.bed

        # Get distance to CDR for every SV call
        bedtools closest \
          -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_snvs.w${WINDOWSIZE}.srt.bed \
          -b /private/groups/migalab/jmmenend/HPRC/cenSatProject/CDR_data/${chr}/${chr}.centrodip.bed \
          -D a -t first \
          > /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/SNVs_pairwise_raw/${WINDOWSIZE}/${chr}/${REF_SMP}_${QUERY_SMP}.CDR_dist.${WINDOWSIZE}.bed
    done

    rm -rf ${LOCAL_FOLDER}
```



### Run for SNVs - filtered csv files

Convert Julian's SNV csv files into bed files in asm coordinates

```sh
#!/bin/bash
#SBATCH --job-name=convert_SNV_bed_to_asm_coords_filt
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --array=[0-23]%24
#SBATCH --time=12:00:00

source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate ipynb

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

chr=${chromosomes[$SLURM_ARRAY_TASK_ID]}

mkdir -p /private/groups/patenlab/mira/centrolign/analysis/SNVs_pairwise_asm_coords_filt/${chr}/

time python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/convert_SNV_csv_to_bed.py \
  -i /private/groups/migalab/juklucas/centrolign/variant_calling/rates/${chr}/snv_calls/filtered/ \
  -a /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/ \
  -c ${chr} \
  -o /private/groups/patenlab/mira/centrolign/analysis/SNVs_pairwise_asm_coords_filt/${chr}/
```

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

chr=${chromosomes[$SLURM_ARRAY_TASK_ID]}

WINDOWSIZE="1kb"

LOCAL_FOLDER=/data/tmp/$(whoami)/${chr}_SNVs_tmps/
mkdir -p ${LOCAL_FOLDER}
mkdir -p /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/SNVs_pairwise_filt/${WINDOWSIZE}/${chr}

cd /private/groups/patenlab/mira/centrolign/analysis/SNVs_pairwise_asm_coords_filt/${chr}/

ls *.bed | while read -r bed ; do
        # for each bed get ref and query sample ID
        REF_SMP=`echo $bed | basename $bed | cut -f1 -d"_"`
        QUERY_SMP=`echo $bed | basename $bed | cut -f2 -d"_" | cut -f1-2 -d"."`

        # get asm contig name
        SMP1_contig=`head -n1 $bed | cut -f1`
        SMP2_contig=`tail -n1 $bed | cut -f1`

        # get file identifier - remove extension
        FILE_ID=${bed%.*}

        # concatenate window files for both samples, restricting just to contigs for this chrom
        grep $SMP1_contig /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/per_smp_asat_beds_windows/${WINDOWSIZE}/${REF_SMP}_asat_arrays.bed > ${LOCAL_FOLDER}/${FILE_ID}.w${WINDOWSIZE}.bed

        grep $SMP2_contig /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/per_smp_asat_beds_windows/${WINDOWSIZE}/${QUERY_SMP}_asat_arrays.bed >> ${LOCAL_FOLDER}/${FILE_ID}.w${WINDOWSIZE}.bed

        # Use bedtools coverage to calculate number of SNV bases overlapping windows
        bedtools coverage \
          -a ${LOCAL_FOLDER}/${FILE_ID}.w${WINDOWSIZE}.bed \
          -b ${bed} \
          > ${LOCAL_FOLDER}/${FILE_ID}.all_snvs.w${WINDOWSIZE}.bed

        # use bedtools closest to get distance from CDR
        bedtools sort -i ${LOCAL_FOLDER}/${FILE_ID}.all_snvs.w${WINDOWSIZE}.bed > ${LOCAL_FOLDER}/${FILE_ID}.all_snvs.w${WINDOWSIZE}.srt.bed

        # Get distance to CDR for every SV call
        bedtools closest \
          -a ${LOCAL_FOLDER}/${FILE_ID}.all_snvs.w${WINDOWSIZE}.srt.bed \
          -b /private/groups/migalab/jmmenend/HPRC/cenSatProject/CDR_data/${chr}/${chr}.centrodip.bed \
          -D a -t first \
          > /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/SNVs_pairwise_filt/${WINDOWSIZE}/${chr}/${FILE_ID}.CDR_dist.${WINDOWSIZE}.bed
    done

    rm -rf ${LOCAL_FOLDER}
```

### Need to show SNV and short indel rate over aligned bases. For every window, per pairwise comparison, get # of aligned bases in the window

Get list of sample pairs per chromosome < 0.2 distance, with their contig IDs to pass to the script
```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices

ls | while read line ; do
    awk -F',' '$3 <= 0.2' $line > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices_lt0.2/$line
done

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

### need to add cigar string path location to contig maps

for chr in "${chromosomes[@]}";
    do
      cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices_lt0.2/${chr}_r2_QC_v2_centrolign_pairwise_distance.csv | while IFS=',' read -r smp1 smp2 dist ; do
        smp1_contig=`grep "${chr}$" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/${smp1}_asat_arrays.bed | cut -f1`
        smp2_contig=`grep "${chr}$" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/${smp2}_asat_arrays.bed | cut -f1`

        echo ${smp1},${smp1_contig},${smp2},${smp2_contig} >> /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/aligned_bases_per_bed/dist_0.2_smp_contig_maps/${chr}.contig_maps.csv
    done
done
```

Run script per chromosome
```sh
#!/bin/bash
#SBATCH --job-name=CDR_dist_snv_windows_raw
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

chr=${chromosomes[$SLURM_ARRAY_TASK_ID]}

mkdir -p /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/aligned_bases_per_bed/SNVs_pairwise_raw_1kb/${chr}

python /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/get_aligned_bases_bed.py \
    -c /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/aligned_bases_per_bed/dist_0.2_smp_contig_maps_cigars/${chr}.contig_maps.csv \
    -b /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/SNVs_pairwise_raw/1kb/${chr}/ \
    -s CDR_dist.1kb \
    -o /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/aligned_bases_per_bed/SNVs_pairwise_raw_1kb/${chr}
```

python /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/get_aligned_bases_bed.py \
    -c /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/aligned_bases_per_bed/test/contig_maps.csv \
    -b /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/aligned_bases_per_bed/test/ \
    -s CDR_dist.1kb \
    -o /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/aligned_bases_per_bed/test/output/
