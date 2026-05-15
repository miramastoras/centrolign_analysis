### Code for processing genes 

1. Get bed file of censat regions for each sample (check Fig1_gene_analysis.ipynb)

2. Subset gff3 to censat regions, dropping ENS genes.
```sh
mkdir -p logs

sbatch /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/filter_gff3.sh /private/groups/patenlab/mira/centrolign/analysis/Fig1_genes/extract_censat_from_gff3/all_chroms_out/
```

3. Annotate GFF3 with censat repeat categories 
```sh
RESULTS=/private/groups/patenlab/mira/centrolign/analysis/Fig1_genes/extract_censat_from_gff3/all_chroms_out/results
CENSAT_DIR=/private/groups/migalab/juklucas/censat_regions/censat_beds
OUTDIR=/private/groups/patenlab/mira/centrolign/analysis/Fig1_genes/annotate_gff3_with_censat/
mkdir -p ${OUTDIR} ${OUTDIR}/results ${OUTDIR}/logs

echo -e "sample_id\thaplotype\tgff3\tcensat_bed" > ${OUTDIR}/jobs.tsv

for gff3 in ${RESULTS}/*.genes.gff3; do
    basename=$(basename $gff3)
    sample_id=$(echo $basename | cut -d'.' -f1)
    haplotype=$(echo $basename | cut -d'.' -f2)

    if [ "$haplotype" == "1" ]; then alt_label="pat"; else alt_label="mat"; fi

    found=$(ls ${CENSAT_DIR}/${sample_id}_*.cenSat.bed 2>/dev/null | grep -E "_hap${haplotype}_|_${alt_label}_" | head -1)

    if [ -z "$found" ]; then
        echo "MISSING: ${sample_id} hap${haplotype}"
        continue
    fi

    echo -e "${sample_id}\t${haplotype}\t${gff3}\t${found}" >> ${OUTDIR}/jobs.tsv
done

wc -l ${OUTDIR}/jobs.tsv
```

Submit slurm script 
```sh
OUTDIR=/private/groups/patenlab/mira/centrolign/analysis/Fig1_genes/annotate_gff3_with_censat/
cd ${OUTDIR}
mkdir -p logs
sbatch /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/annotate_censat_genes.sh ${OUTDIR}
```

Fix gff3 for IGV 
awk 'BEGIN{OFS="\t"} /^#/{print; next} {$1="HG01071#1#"$1; print}' HG01071.1.genes.gff3 > HG01071.1.genes.fix.gff3
awk 'BEGIN{OFS="\t"} /^#/{print; next} {$1="HG00408#1#"$1; print}' HG00408.1.genes.gff3 > /private/groups/patenlab/mira/HG00408.1.genes.fix.gff3
awk 'BEGIN{OFS="\t"} /^#/{print; next} {$1="HG00097#2#"$1; print}' HG00097.2.genes.gff3 > /private/groups/patenlab/mira/HG00097.2.genes.fix.gff3


### starting over with Prajna's gene list 

```
/private/groups/patenlab/mira/centrolign/analysis/Fig1_genes/annotate_censat_prajnas_list/

bash /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/annotate_censat_genes_intersect.sh --make-jobs

sbatch /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/annotate_censat_genes_intersect.sh
```