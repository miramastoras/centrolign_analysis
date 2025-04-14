## Preparing list of assemblies

Identify samples that align well to other samples in Faith's test set GFA, but are not in the GFA

```sh
# all pairwise alignments in release 2 with < 0.5 alignment distance
awk -F, '$3 < 0.5' /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr12/pairwise_distance.csv  > /private/groups/patenlab/mira/centrolign/giraffe/pairwise_dist_lt0.5.csv

# all pairwise alignments involving samples inside Faith's test graph
cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/fasta_list.all_sample_ids.txt | while read line ; do grep $line /private/groups/patenlab/mira/centrolign/giraffe/pairwise_dist_lt0.5.csv ; done | sort | uniq > /private/groups/patenlab/mira/centrolign/giraffe/pairwise_alns_inside_graph.csv

# get all sample names involved in these pairwise alignments
cut -f 1 -d"," /private/groups/patenlab/mira/centrolign/giraffe/pairwise_alns_inside_graph.csv > samps
cut -f 2 -d"," /private/groups/patenlab/mira/centrolign/giraffe/pairwise_alns_inside_graph.csv >> samps

sort samps | uniq > tmp ; mv tmp samps

# now sample names that are NOT in graph
grep -v -f /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/fasta_list.all_sample_ids.txt samps > /private/groups/patenlab/mira/centrolign/giraffe/samples_aligning_well_not_in_graph.txt

# get their fasta locations
cat /private/groups/patenlab/mira/centrolign/giraffe/samples_aligning_well_not_in_graph.txt | while read line ; do grep $line /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/HPRC_release2_contiguous_HORs_chr12.fasta_list.txt ; done > /private/groups/patenlab/mira/centrolign/giraffe/fastas_not_in_graph_to_align.txt

# append full path
sed 's/^/\/private\/groups\/patenlab\/mira\/centrolign\/batch_submissions\/extract_hors_HPRC\/release2\//' /private/groups/patenlab/mira/centrolign/giraffe/fastas_not_in_graph_to_align.txt > tmp ; mv tmp /private/groups/patenlab/mira/centrolign/giraffe/fastas_not_in_graph_to_align.txt
```
## Extracting HiFi read alignments

Sample in Faith's graph: `HG00099.1`
chosen from here:
```
/private/groups/patenlab/mira/centrolign/giraffe/pairwise_alns_inside_graph.csv
```

Sample not in Faith's graph which aligns well to others in the graph `HG00706.1`
chosen from here:

```
/private/groups/patenlab/mira/centrolign/giraffe/fastas_not_in_graph_to_align.txt
```

Data table with locations for the alignments of HiFi reads to the assemblies [located here](https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/data_tables/assembly_qc/flagger/flagger_hifi_processing_metadata_v0.1.csv)


Download HiFi reads for HG00706 and HG00099 to phoenix:
```sh
aws s3 cp s3://human-pangenomics/submissions/ca366a13-5bad-487b-8a57-97344e9aa0e4--HPRC_RELEASE_2_SUPPLEMENTARY_ASSEMBLY_QC/HG00099/hprc_r2/assembly_qc/read_alignments/hifi/HG00099.hifi_DC_minimap2_2.28.corrected.bam /private/groups/patenlab/mira/centrolign/giraffe/extract_hifi_reads/

aws s3 cp s3://human-pangenomics/submissions/ca366a13-5bad-487b-8a57-97344e9aa0e4--HPRC_RELEASE_2_SUPPLEMENTARY_ASSEMBLY_QC/HG00706/hprc_r2/assembly_qc/read_alignments/hifi/HG00706.hifi_minimap2_2.28.corrected.bam /private/groups/patenlab/mira/centrolign/giraffe/extract_hifi_reads/
```
Get coordinates of the chr12 HOR array for HG00099.1 and HG00706.1

> .1 means hap1 or pat, .2 means hap2 or mat

```
grep chr12 /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/contiguous_HORs_bed_files/HG00706_pat_hprc_r2_v1.0.1_hor_arrays.bed \
> /private/groups/patenlab/mira/centrolign/giraffe/extract_hifi_reads/HG00706_pat_hprc_r2_v1.0.1_hor_arrays.chr12.bed

grep chr12 /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/contiguous_HORs_bed_files/HG00099_hap1_hprc_r2_v1_hor_arrays.bed \
> /private/groups/patenlab/mira/centrolign/giraffe/extract_hifi_reads/HG00099_hap1_hprc_r2_v1_hor_arrays.chr12.bed
```

Extract reads overlapping the HOR arrays
```sh
#!/bin/bash
#SBATCH --job-name=extract_reads
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=samtools_%x.%j.log
#SBATCH --time=1:00:00

samtools view -bh -@32 \
  -L /private/groups/patenlab/mira/centrolign/giraffe/extract_hifi_reads/HG00706_pat_hprc_r2_v1.0.1_hor_arrays.chr12.bed \
  /private/groups/patenlab/mira/centrolign/giraffe/extract_hifi_reads/HG00706.hifi_minimap2_2.28.corrected.bam \
  > /private/groups/patenlab/mira/centrolign/giraffe/extract_hifi_reads/HG00706_pat_chr12_hor_array.hifi.bam

samtools index /private/groups/patenlab/mira/centrolign/giraffe/extract_hifi_reads/HG00706_pat_chr12_hor_array.hifi.bam

samtools view -bh -@32 \
  -L /private/groups/patenlab/mira/centrolign/giraffe/extract_hifi_reads/HG00099_hap1_hprc_r2_v1_hor_arrays.chr12.bed \
  /private/groups/patenlab/mira/centrolign/giraffe/extract_hifi_reads/HG00099.hifi_DC_minimap2_2.28.corrected.bam \
  > /private/groups/patenlab/mira/centrolign/giraffe/extract_hifi_reads/HG00099_hap1_chr12_hor_array.hifi.bam

samtools index /private/groups/patenlab/mira/centrolign/giraffe/extract_hifi_reads/HG00099_hap1_chr12_hor_array.hifi.bam
```
