## Submitting all vs all centrolign direct pairwise alignments for release2 QC v2

"release2 QCv2" indicates samples passing Julian's QC pipeline, which integrates Flagger and
nucFlag to filter arrays.


Downloaded full assembly list for release 2 from here https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/data_tables/assemblies_pre_release_v0.6.1.index.csv

```sh
/private/groups/patenlab/mira/centrolign/annotations/assemblies_pre_release_v0.6.1.index.csv
```

### Starting with Chr 1, Chr 12, Chr 10, Chr 11, Chr 8


#### 1. Extract HOR fasta files

```sh
# Add s3 location to csv and get table in correct format
chromosomes=("chr1" "chr12" "chr10" "chr11" "chr8")

for chr in "${chromosomes[@]}"
do
    mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/active_array_fastas/${chr}/

    cat /private/groups/migalab/juklucas/censat_regions/active_arrays/asat_arrays_${chr}.csv | grep -v "sample_id" | while read line ; do
        asm_name=`echo $line | cut -f3 -d","`
        fasta_links=`grep $asm_name /private/groups/patenlab/mira/centrolign/annotations/assemblies_pre_release_v0.6.1.index.csv | cut -f11,13 -d","`

        echo "$line" | awk -F, -v fasta_links="$fasta_links" 'BEGIN {OFS=","} {print $1, $2, $3, $4, $5, $6, $16, fasta_links}' >> /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/active_array_csvs/${chr}_asat_arrays_r2_QCv2.csv
    done
done

### slurm script to extract fasta

# extract fasta files
samtools faidx -r $REGIONFILE ~{assemblyFasta} | sed "s/>/>$SAMPLE.$PARNUM /g" > $HORFASTA
```

```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/parse_QC_csv.py /private/groups/migalab/juklucas/censat_regions/active_arrays/asat_arrays_chr1.csv /private/groups/patenlab/mira/centrolign/annotations/assemblies_pre_release_v0.6.1.index.csv
```
