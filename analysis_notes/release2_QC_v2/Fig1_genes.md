### Code for processing genes 

1. Get bed file of censat regions for each sample (check Fig1_gene_analysis.ipynb)

2. Subset gff3 to censat regions, dropping ENS genes.
```sh
mkdir -p logs

sbatch /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/filter_gff3.sh /private/groups/patenlab/mira/centrolign/analysis/Fig1_genes/extract_censat_from_gff3/all_chroms_out/
```