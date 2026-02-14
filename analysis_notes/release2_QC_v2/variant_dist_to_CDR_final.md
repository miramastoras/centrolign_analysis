## Calculating CDR distance to mutation windows for updated centrodip calls (2/10/26)

### SNVs


SNVs for 10kb:
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

WINDOWSIZE=$1

LOCAL_FOLDER=/data/tmp/$(whoami)/${chr}_SNVs_${WINDOWSIZE}_tmps/
mkdir -p ${LOCAL_FOLDER}

OUTDIR=/private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR_centrodip_02102026/SNVs_pairwise_raw/${WINDOWSIZE}/${chr}
mkdir -p $OUTDIR

cd /private/groups/patenlab/mira/centrolign/analysis/SNVs_induced_pairwise_asm_coords/${chr}/

cat /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/bed_lists_dist_0.2/SNVs_induced_raw/${chr}.bed_paths.txt | while read -r bed ; do
        # for each bed get ref and query sample ID
        REF_SMP=`echo $bed | basename $bed | cut -f1 -d"_"`
        QUERY_SMP=`echo $bed | basename $bed | cut -f2 -d"_" | cut -f1-2 -d"."`

        # convert bedPE into bed file, containing SV calls for both samples
        cut -f1-3 $bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_SNVs.bed
        cut -f4-6 $bed >> ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_SNVs.bed

        # get asm contig name
        SMP1_contig=`head -n1 ${bed} | cut -f1`
        SMP2_contig=`head -n1 ${bed} | cut -f4`

        # concatenate window files for both samples, restricting just to contigs for this chrom
        grep $SMP1_contig /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/per_smp_asat_beds_windows/${WINDOWSIZE}/${REF_SMP}_asat_arrays.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.bed

        grep $SMP2_contig /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/per_smp_asat_beds_windows/${WINDOWSIZE}/${QUERY_SMP}_asat_arrays.bed >> ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.bed

        # Use bedtools coverage to calculate number of SNV bases overlapping windows
        bedtools coverage \
          -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.w${WINDOWSIZE}.bed \
          -b ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_SNVs.bed \
          > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_snvs.w${WINDOWSIZE}.bed

        # use bedtools closest to get distance from CDR
        bedtools sort -i ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_snvs.w${WINDOWSIZE}.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_snvs.w${WINDOWSIZE}.srt.bed

        # Get distance to CDR for every SV call
        bedtools closest \
          -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_snvs.w${WINDOWSIZE}.srt.bed \
          -b /private/groups/patenlab/mira/centrolign/analysis/CDR_permutation/per_chrom_centrodip_2_10_26/${chr}_CDRs.bed.sort.bed \
          -D a -t first \
          > ${OUTDIR}/${REF_SMP}_${QUERY_SMP}.CDR_dist.${WINDOWSIZE}.bed
    done

rm -rf ${LOCAL_FOLDER}
```

```sh
cd /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR_centrodip_02102026/SNVs_pairwise_raw/

mkdir -p logs
sbatch variant_dist_CDR.sh "100kb"
sbatch variant_dist_CDR.sh "10kb"
```

### short indels

```sh
#!/bin/bash
#SBATCH --job-name=CDR_dist_short_indel_windows
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/call_SVs_%x.%j.log
#SBATCH --time=1:00:00
#SBATCH --array=[0-23]%24

# activate environment for bedtools
source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate base

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

chr=${chromosomes[$SLURM_ARRAY_TASK_ID]}

WINDOWSIZE=$1

LOCAL_FOLDER=/data/tmp/$(whoami)/${chr}_short_indels_tmp_${WINDOWSIZE}/
mkdir -p ${LOCAL_FOLDER}
OUTDIR=/private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR_centrodip_02102026/short_indels_pairwise/${WINDOWSIZE}/${chr}
mkdir -p ${OUTDIR}

cat /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/bed_lists_dist_0.2/short_indels_induced/${chr}.bed_paths.txt | while read -r bed ; do

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
          -b /private/groups/patenlab/mira/centrolign/analysis/CDR_permutation/per_chrom_centrodip_2_10_26/${chr}_CDRs.bed.sort.bed \
          -D a -t first \
          > ${OUTDIR}/${REF_SMP}_${QUERY_SMP}.CDR_dist.${WINDOWSIZE}.bed
    done

rm -rf ${LOCAL_FOLDER}

```

```sh
cd /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR_centrodip_02102026/short_indels_pairwise

mkdir -p logs
sbatch CDR_dist.sh "10kb"
sbatch CDR_dist.sh "100kb"
```


### SVs

```sh
#!/bin/bash
#SBATCH --job-name=CDR_dist_SVs_windows
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x.%j.log
#SBATCH --time=1:00:00
#SBATCH --array=[3]%24

# activate environment for bedtools
source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate base

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

chr=${chromosomes[$SLURM_ARRAY_TASK_ID]}

WINDOWSIZE=$1

LOCAL_FOLDER=/data/tmp/$(whoami)/${chr}_SVs_${WINDOWSIZE}_tmp/
mkdir -p ${LOCAL_FOLDER}
OUTDIR=/private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR_centrodip_02102026/SVs_pairwise/${WINDOWSIZE}/${chr}
mkdir -p ${OUTDIR}

cat /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/bed_lists_dist_0.2/SVs_induced/${chr}.bed_paths.txt | while read -r bed ; do

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
          -b /private/groups/patenlab/mira/centrolign/analysis/CDR_permutation/per_chrom_centrodip_2_10_26/${chr}_CDRs.bed.sort.bed \
          -D a -t first \
          > ${OUTDIR}/${REF_SMP}_${QUERY_SMP}.CDR_dist.${WINDOWSIZE}.bed
    done

rm -rf ${LOCAL_FOLDER}
```

Submitting
```sh
cd /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR_centrodip_02102026/SVs_pairwise/

mkdir -p logs
sbatch CDR_dist_SVs_windows.sh "100kb"
sbatch CDR_dist_SVs_windows.sh "50kb"
```


#### Add aligned bases per window into the bed file

Run script per chromosome
```sh
#!/bin/bash
#SBATCH --job-name=get_aligned_bases_short_indels
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

chr=${chromosomes[$SLURM_ARRAY_TASK_ID]}

WINDOWSIZE=$1
OUTDIR=/private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR_centrodip_02102026/aligned_bases_per_bed/short_indels_pairwise_${WINDOWSIZE}/${chr}
mkdir -p $OUTDIR

python /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/get_aligned_bases_bed.py \
    -c /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/aligned_bases_per_bed/dist_0.2_smp_contig_maps_cigars/${chr}.contig_maps.csv \
    -b /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR_centrodip_02102026/short_indels_pairwise/${WINDOWSIZE}/${chr}/ \
    -s CDR_dist.${WINDOWSIZE} \
    -o ${OUTDIR}
```

```sh
cd /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR_centrodip_02102026/aligned_bases_per_bed

mkdir -p logs
sbatch short_indels.sh "10kb"
sbatch short_indels.sh "100kb"
```


SNVs

```sh
#!/bin/bash
#SBATCH --job-name=get_aligned_bases_snvs
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

WINDOWSIZE=$1
OUTDIR=/private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR_centrodip_02102026/aligned_bases_per_bed/SNVs_pairwise_${WINDOWSIZE}/${chr}
mkdir -p $OUTDIR

python /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/get_aligned_bases_bed.py \
    -c /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR/aligned_bases_per_bed/dist_0.2_smp_contig_maps_cigars/${chr}.contig_maps.csv \
    -b /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR_centrodip_02102026/SNVs_pairwise_raw/${WINDOWSIZE}/${chr}/ \
    -s CDR_dist.${WINDOWSIZE} \
    -o ${OUTDIR}
```


```sh
cd /private/groups/patenlab/mira/centrolign/analysis/variant_dist_CDR_centrodip_02102026/aligned_bases_per_bed

mkdir -p logs
sbatch snvs.sh "10kb"
sbatch snvs.sh "50kb"
sbatch snvs.sh "100kb"
```
