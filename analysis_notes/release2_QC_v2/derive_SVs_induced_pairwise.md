## Deriving SV calls from induced pairwise cigar strings

### Starting with Chr8 as a test case

#### Step 1: Identify low divergence clades to call variants in

> Using Julian's default of 95% of samples having < 0.8 distance

chr8 subgroup 1 - all one clade

```sh
docker run -u `id -u`:`id -g` -v /private/groups:/private/groups/ \
    miramastoras/centromere_scripts:v0.1.4 \
    python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/identify_low_divergence_clades.py \
  --newick /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr8/split_by_2/subgroup_1_tree.nwk \
  --distance_csv /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/chr8_r2_QC_v2_centrolign_pairwise_distance.csv \
  --max_pairwise_dist 0.8 \
  --output_prefix /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/subgroup_1
```
chr 8 subgroup 0 - 8 clades
```sh
docker run -u `id -u`:`id -g` -v /private/groups:/private/groups/ \
    miramastoras/centromere_scripts:v0.1.4 \
    python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/identify_low_divergence_clades.py \
  --newick /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr8/split_by_2/subgroup_0_tree.nwk \
  --distance_csv /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/chr8_r2_QC_v2_centrolign_pairwise_distance.csv \
  --max_pairwise_dist 0.8 \
  --output_prefix /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/subgroup_0
```

Color subgroup 0 clade 2 (largest one)
```sh
chromosomes=("chr8")

for chr in "${chromosomes[@]}"
do
  python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap_v2.py \
    -t /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_all_pairs_nj_tree.format5.nwk \
    -s /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}.samples.txt \
    -p /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_pairwise_distance.csv \
    -m "Centrolign all pairs distances" \
    -n "${chr} NJ tree" \
    -d "All pairs Distances" \
    -o /Users/miramastoras/Desktop/color_subgroups_heatmap/${chr}_r2_QC_v2_clade_2_subgroup0 --no_labels \
    --highlight_samples /Users/miramastoras/Desktop/color_subgroups_heatmap/chr8_subgroup0_clade2.txt
done
```

#### Step 2: Call SVs from cigar strings within each clade

input: directory containing cigar strings
output: bedPE file with SV coords relative to each assembly

Testing how long it takes to run on chr8 subgroup 1
```sh
cut -f1 subgroup_1_seqs.fasta.fai  > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/chr8_subgroup1_samples.txt
```

```sh
#!/bin/bash
#SBATCH --job-name=SVs_chr8_subgroup_1
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=call_SVs_%x.%j.log
#SBATCH --time=7-00:00

time python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/call_SVs_pairwise.py \
  -c /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr8/subgroup_1/induced_pairwise_cigars/pairwise_cigar_ \
  -s /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/chr8_subgroup1_samples.txt \
  -o /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/subgroup_1_SV_beds/
```
Ran in two minutes

Local computer testing:
```sh
time python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/call_SVs_pairwise.py \
  -c /Users/miramastoras/Desktop/chr12_cigars/pairwise_cigar_ \
  -s /Users/miramastoras/Desktop/chr12_cigars/samples_test.txt \
  -o /Users/miramastoras/Desktop/chr12_cigars/SVs_
```
Plot SV length distributions
```sh
docker run -u `id -u`:`id -g` -v /private/groups:/private/groups/ \
    miramastoras/centromere_scripts:v0.1.4 \
    python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/plot_SVs_pairwise.py \
    -b /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/subgroup_1_SV_beds/ \
    -o /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/SV_violin_subgroup1.png
```
Total SVs: 2409434 ({'D': 1204717, 'I': 1204717})

#### Run SV caller on all chr 8 subgroups

Create directory containing sample lists for all chr 8 subgroups
```sh
# subgroup 1
cut -f1 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr8/split_by_2/subgroup_1_seqs.fasta.fai  > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/clade_sample_lists/subgroup1_samples.txt

# subgroup 0
clades=("Clade_1" "Clade_2" "Clade_3" "Clade_4" "Clade_5" "Clade_6" "Clade_7" "Clade_8")

for clade in "${clades[@]}"
do
  grep $clade /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/subgroup_0_clades.csv | cut -f2 -d"," > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/clade_sample_lists/subgroup_0_${clade}_samples.txt
done
```
Get list of clade names to run:
```sh
cd /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/clade_sample_lists

for f in *_samples.txt; do
  echo "${f%_samples.txt}"
done > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/clade_names.txt

# manually added cigar string path to column 2
```

Run SV caller on all subgroups
```sh
#!/bin/bash
#SBATCH --job-name=SVs_chr8
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/call_SVs_%x.%j.log
#SBATCH --time=1:00:00
#SBATCH --array=[1-9]%9

CLADE_NAMES=/private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/clade_names.txt
CLADE=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$CLADE_NAMES" | cut -f1 -d",")
CIGAR_PATH=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$CLADE_NAMES" | cut -f2 -d",")

echo $CLADE
echo $CIGAR_PATH

mkdir -p /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/SV_beds/${CLADE}/

time python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/call_SVs_pairwise.py \
  -c ${CIGAR_PATH} \
  -s /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/clade_sample_lists/${CLADE}_samples.txt \
  -o /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/SV_beds/${CLADE}/
```

```sh
cd /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/

mkdir -p logs
mkdir -p /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/SV_beds/
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/ pull
sbatch slurm_call_SVs.sh
```

Create csv file for plotting results
```sh
cd /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test

echo "clade,path_to_SV_beds" > plotting_script_samples.csv
cat /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/clade_names.txt | cut -f1 -d"," | while read line ; do
  echo $line,"/private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/SV_beds/${line}"
done >> plotting_script_samples.csv
```

Generate plots for chr8
```sh
time python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/plot_SVs_pairwise.py \
  -i /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/plotting_script_samples.csv \
  -p /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr8/test/plots/
```
Transitioned to python notebook for remaining plots:
https://github.com/miramastoras/centrolign_analysis/blob/main/analysis_notes/release2_QC_v2/notebooks/SVs_pairwise.ipynb

### Running SV caller on all pairwise comparisons for all chromosomes

Get lists containing all the samples in each subgroup / MSA submission
```sh
grep "complete" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/centrolign_MSA.csv | grep "release 2 QC v2" | cut -f1,8 -d","  | while IFS=',' read -r subgroup fastapath ; do
  CHR=`echo $subgroup | cut -f1 -d"_"`
  echo $CHR
  mkdir -p /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/${CHR}/

  cut -f1 ${fastapath}.fai > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/sample_lists/${subgroup}.MSA.samples.txt
done
```
Create csv file with columns clade,CIGAR_PATH,sample_names_list
```sh
grep "complete" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/centrolign_MSA.csv | grep "release 2 QC v2" | cut -f1,8,10 -d","  > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/11172025_completed_subgroups.csv
```
Run SV caller on all completed subgroups
```sh

```

submitting
```sh
cd /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise

mkdir -p logs

sbatch /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/call_SVs_pairwise_all_chroms.sh /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/11172025_completed_subgroups.csv
```
### Benchmarking: Compare Centrolign to Fedor's HorHap SVs

Starting with chr 12 cenHap 4.

```sh
/private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps
```
On local computer, plot cenhap 4 with respect to chr12 HOR tree - confirmed they are in subgroup 1
```sh
chromosomes=("chr12")

for chr in "${chromosomes[@]}"
do
  python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap_v2.py \
    -t /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_all_pairs_nj_tree.format5.nwk \
    -s /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}.samples.txt \
    -p /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_pairwise_distance.csv \
    -m "Centrolign all pairs distances" \
    -n "${chr} NJ tree" \
    -d "All pairs Distances" \
    -o /Users/miramastoras/Desktop/color_subgroups_heatmap/${chr}_r2_QC_v2_cenHap4_colored --no_labels \
    --highlight_samples /Users/miramastoras/Desktop/chr12_cenhap4_HPRC_sample_list.txt
done
```
Run SV calling script on cenhap 4 samples
```sh
#!/bin/bash
#SBATCH --job-name=SVs_chr12_cenhap4
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/call_SVs_%x.%j.log
#SBATCH --time=1:00:00


mkdir -p /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/SV_beds/

time python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/call_SVs_pairwise.py \
  -c /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/induced_pairwise_cigars/pairwise_cigar_ \
  -s /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/chr12_cenhap4_HPRC_sample_list.txt \
  -o /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/SV_beds/
```

Convert coordinates back to asm coords
```sh
time python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/convert_SV_bed_to_asm_coords.py \
  -s /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/SV_beds/ \
  -a /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/ \
  -c chr12 \
  -o /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/SV_beds_asm_coords/
```

Comparing percent overlap with Fedor's SVs.


Ran this python notebook to convert bed files into comparable format: https://github.com/miramastoras/centrolign_analysis/blob/main/analysis_notes/release2_QC_v2/notebooks/SVs_pairwise.ipynb

```sh
# report all centrolign SVs and their percent overlap with horhap SVs
bedtools intersect \
    -a /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/centrolign_SVs_ins.bed \
    -b /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/horhap_SVs_ins.bed \
    -wao | awk '{pct = ($NF / ($3 - $2)) * 100; print $0 "\t" pct}'


# Calculate percent of centrolign features overlapped by HorHap features - ins

bedtools sort -i /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/horhap_SVs_ins.bed | bedtools merge -i - > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/horhap_SVs_ins.merged.bed

bedtools coverage \
    -a /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/centrolign_SVs_ins.bed \
    -b /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/horhap_SVs_ins.merged.bed \
    > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/bedtools_coverage_ins.bed

# Calculate percent of centrolign features overlapped by HorHap features - del

bedtools sort -i /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/horhap_SVs_del.bed | bedtools merge -i - > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/horhap_SVs_del.merged.bed

bedtools coverage \
    -a /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/centrolign_SVs_del.bed \
    -b /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/horhap_SVs_del.merged.bed \
    > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/bedtools_coverage_del.bed

# combine bed file
cat /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/bedtools_coverage_del.bed  /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/bedtools_coverage_ins.bed >  /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/bedtools_coverage_del_ins.bed


# Calculate percent of HorHap features overlapped by Centrolign features - ins

bedtools sort -i /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/centrolign_SVs_ins.bed | bedtools merge -i - > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/centrolign_SVs_ins.mrg.bed

bedtools coverage \
    -b /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/centrolign_SVs_ins.mrg.bed \
    -a /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/horhap_SVs_ins.bed \
    > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/bedtools_coverage_ins_rev.bed

# Calculate percent of centrolign features overlapped by HorHap features - del

bedtools sort -i /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/centrolign_SVs_del.bed | bedtools merge -i - > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/centrolign_SVs_del.mrg.bed


bedtools coverage \
    -b /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/centrolign_SVs_del.mrg.bed \
    -a /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/horhap_SVs_del.bed \
    > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/bedtools_coverage_del_rev.bed

# combine bed file
cat /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/bedtools_coverage_del_rev.bed  /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/bedtools_coverage_ins_rev.bed >  /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/bedtools_coverage_del_ins_rev.bed
```

Using the tool "intervene" to plot a venn diagram of bedtools intersect
```sh
conda create -n intervene
conda install -c bioconda intervene

intervene venn -i /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/centrolign_SVs_ins.bed /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/horhap_SVs_ins.bed –bedtools-options -f 0.5
```
