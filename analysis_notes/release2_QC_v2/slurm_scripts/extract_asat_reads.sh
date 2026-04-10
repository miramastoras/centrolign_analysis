#!/bin/bash
#SBATCH --job-name=extract_ASAT_reads
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=56gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --array=[2-232]%50
#SBATCH --exclude=phoenix-[09,10,22,23,24,18]
#SBATCH --output=/private/groups/patenlab/mira/centrolign/rGFA_tests/extract_asat_reads/logs/array_job_%A_task_%a.log
#SBATCH --time=1:00:00

READ_LOCS=/private/groups/patenlab/fokamoto/centrolign/to_align/aws_file_locations.csv
BED_DIR=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds
OUT_BASE=/private/groups/patenlab/mira/centrolign/rGFA_tests/extract_asat_reads/

SAMPLE_ID=`head -n "$SLURM_ARRAY_TASK_ID" "$READ_LOCS" | tail -n 1 | cut -f1 -d ","`
echo "Running sample: $SAMPLE_ID"

LOCAL_FOLDER=/data/tmp/$(whoami)/${SAMPLE_ID}_asat_extract/
mkdir -p ${LOCAL_FOLDER}

echo "Downloading reads for $SAMPLE_ID from AWS"
reads=`grep "^$SAMPLE_ID," "$READ_LOCS" | cut -f3 -d ","`
reads_index=$(grep "^$SAMPLE_ID," "$READ_LOCS" | cut -f4 -d ",")
aws s3 --no-sign-request cp "$reads" $LOCAL_FOLDER/full.bam &> /dev/null
aws s3 --no-sign-request cp "$reads_index" ${LOCAL_FOLDER}/full.bam.bai &> /dev/null
echo "Download complete"

for hap_num in 1 2; do
    hap_id="${SAMPLE_ID}.${hap_num}"

    # Subset all-chr centromere coordinate file to have only our chosen ones
    # -e chr4 -e chr6 -e chr9 -e chr10 -e chr11 -e chr12 -e chr17 --> already extracted by Faith
    grep -e 'chr1$' -e 'chr2$' -e 'chr3$' -e 'chr5$' -e 'chr7$' -e 'chr8$' \
         -e 'chr13$' -e 'chr14$' -e 'chr15$' -e 'chr16$' -e 'chr18$' -e 'chr19$' \
         -e 'chr20$' -e 'chr21$' -e 'chr22$' -e 'chrX$' -e 'chrY$' \
        ${BED_DIR}/${hap_id}_asat_arrays.bed > ${LOCAL_FOLDER}/chosen_chroms.bed

    if [ -s ${LOCAL_FOLDER}/chosen_chroms.bed ]; then
        echo "BED coordinates found for $hap_id"
        cat ${LOCAL_FOLDER}/chosen_chroms.bed
    else
        echo "No BED coordinates found for $hap_id"
        continue
    fi

    # Subset giant bam to all asats
    samtools view -@32 -L ${LOCAL_FOLDER}/chosen_chroms.bed -b ${LOCAL_FOLDER}/full.bam > ${LOCAL_FOLDER}/centromeres.bam
    samtools index ${LOCAL_FOLDER}/centromeres.bam

    cat ${LOCAL_FOLDER}/chosen_chroms.bed | while read line; do
        echo "$line"
        contig_id=$(echo "$line" | cut -f1)
        start=$(echo "$line" | cut -f2)
        end=$(echo "$line" | cut -f3)
        chrom=$(echo "$line" | cut -f4)

        mkdir -p ${OUT_BASE}/${chrom}/

        echo "Extracting $hap_id reads for $chrom ($contig_id:$start-$end)"

        # Extract reads overlapping this region, rename header contig, convert to fastq.gz
        samtools view -@32 -h ${LOCAL_FOLDER}/centromeres.bam "${contig_id}:${start}-${end}" \
        | samtools fastq -@32 - \
        | gzip > ${OUT_BASE}/${chrom}/${hap_id}_${chrom}_asat.fastq.gz

        if [ -s ${OUT_BASE}/${chrom}/${hap_id}_${chrom}_asat.fastq.gz ]; then
            echo "Created fastq.gz for $hap_id $chrom"
        else
            echo "WARNING: empty fastq.gz for $hap_id $chrom"
        fi
    done

done

rm -rf ${LOCAL_FOLDER}