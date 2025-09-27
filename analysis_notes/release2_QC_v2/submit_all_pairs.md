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
# Python script to create new csv file with s3 links
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

chromosomes=("chr1" "chr12" "chr10" "chr11" "chr8")
for chr in "${chromosomes[@]}"
do
    mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/active_array_fastas/${chr}/

    mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/active_array_beds/${chr}/

    python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/parse_QC_csv.py \
      /private/groups/migalab/juklucas/censat_regions/active_arrays/asat_arrays_${chr}.csv \
      /private/groups/patenlab/mira/centrolign/annotations/assemblies_pre_release_v0.6.1.index.csv \
      /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/active_array_csvs/asat_arrays_${chr}.csv
      /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/active_array_beds/${chr}/
done
```
Run slurm script to extract fasta file
```

```
