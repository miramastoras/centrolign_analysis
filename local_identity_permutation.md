### Permutation analysis for

Get List of local identity beds per sample, converted to asm coordinates

```sh

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

for CHR in "${chromosomes[@]}";
  do

    # index file for local identity beds
    LOCAL_ID_IDX=/private/groups/migalab/juklucas/censat_paper/local_identity/2025_11_13/${CHR}/${CHR}_local_id_launch_outputs.csv

    LOCAL_FOLDER=/private/groups/patenlab/mira/centrolign/analysis/local_identity/permutation/local_id_beds_asm_coords/${CHR}/

    mkdir -p ${LOCAL_FOLDER}

    grep -v "chm13" $LOCAL_ID_IDX | grep -v "hg002v1.1" | cut -f4,11 -d"," | tail -n +2 | while IFS=',' read -r contig local_id_bed ; do

      SAMPLE=`echo $contig | cut -f1-2 -d"#" | sed 's/#/./g'`

      # convert local identity bed to assembly coordinates
      ASAT_START=`grep "${CHR}$" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/${SAMPLE}_asat_arrays.bed | cut -f2`

      awk -F'\t' -v OFS='\t' -v n=$ASAT_START '{ $2 += n; $3 += n }1' ${local_id_bed} > ${LOCAL_FOLDER}/${SAMPLE}_${CHR}_local_id_asat_coords.bed

    done
done
```
