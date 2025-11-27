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
grep "complete" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/centrolign_MSA.csv | grep "release 2 QC v2" | cut -f1,10 -d","  | while IFS=',' read -r subgroup cigar ; do
  echo ${subgroup},${cigar},/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/sample_lists/${subgroup}.MSA.samples.txt >> /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/11172025_completed_subgroups.csv
done
```
Run SV caller on all completed subgroups
```sh
cd /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise

mkdir -p logs

sbatch /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/call_SVs_pairwise_all_chroms.sh /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/11172025_completed_subgroups.csv /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise
```

Create csv file formatted as clade,chr,path_to_SV_beds
```sh
cd /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise

echo "clade,chr,path_to_SV_beds" > 11172025_clade_chr_sv_beds.csv
cat /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/11172025_completed_subgroups.csv | while IFS=',' read -r clade cigars sample_lists ; do
  CHR=$(echo $clade | cut -f1 -d"_")
  echo $clade,$CHR,/private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/${CHR}/SV_beds/${clade}/ >> 11172025_clade_chr_sv_beds.csv
done
```

Get lists of all samples with pairwise distance less than 0.4
```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices

ls | while read line ; do
    CHR=`echo $line | cut -f1 -d"_"`
    echo $CHR
    awk -F',' '$3 <= 0.4' $line | cut -f1-2 -d"," >> /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/pairwise_cutoffs/${CHR}_pairwise_dist.0.4.csv
done
```

Used python notebook to get CSV files from comparisons that would be made from the MSA not direct pairwise


Plot the sample comparisons that would be made
```sh
chromosomes=("chr1" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrY")

for chr in "${chromosomes[@]}" ; do
    { echo "sample1,sample2,dist"; cat "/Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_pairwise_distance.csv"; } > tmp && mv tmp "/Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_pairwise_distance.csv"
    echo "Added header to /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_pairwise_distance.csv"
done

# 0.4
for chr in "${chromosomes[@]}"
do
  python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap_v2.py \
    -t /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_all_pairs_nj_tree.format5.nwk \
    -s /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}.samples.txt \
    -p /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_pairwise_distance.csv \
    -m "Centrolign all pairs distances" \
    -n "${chr} NJ tree" \
    -d "All pairs Distances" \
    -o /Users/miramastoras/Desktop/github_repos/centrolign_analysis/analysis_notes/release2_QC_v2/plots/${chr}.induced_MSA.dist_0.4_ --no_labels \
    --highlight_pairs /Users/miramastoras/Desktop/pairwise_cutoffs/${chr}.induced_MSA.dist_0.4.csv
done

# 0.6
for chr in "${chromosomes[@]}"
do
  python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap_v2.py \
    -t /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_all_pairs_nj_tree.format5.nwk \
    -s /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}.samples.txt \
    -p /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_pairwise_distance.csv \
    -m "Centrolign all pairs distances" \
    -n "${chr} NJ tree" \
    -d "All pairs Distances" \
    -o /Users/miramastoras/Desktop/github_repos/centrolign_analysis/analysis_notes/release2_QC_v2/plots/${chr}.induced_MSA.dist_0.6_ --no_labels \
    --highlight_pairs /Users/miramastoras/Desktop/pairwise_cutoffs/${chr}.induced_MSA.dist_0.6.csv
done
```

### Run SV caller on direct pairwise alignments

Just running this on completed subgroups from MSA for now - using sample lists from above
```sh
/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/sample_lists/${subgroup}.MSA.samples.txt
```
Create csv file with columns clade,CIGAR_PATH,sample_names_list
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/analysis/SVs_direct_pairwise

grep "complete" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/centrolign_MSA.csv | grep "release 2 QC v2" | cut -f1,10 -d","  | while IFS=',' read -r subgroup cigar ; do
  CHR=`echo ${subgroup} | cut -f1 -d"_"`
  echo ${subgroup},/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/${CHR}/pairwise_cigar,/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/sample_lists/${subgroup}.MSA.samples.txt >> /private/groups/patenlab/mira/centrolign/analysis/SVs_direct_pairwise/11172025_completed_subgroups.csv
done
```
Run SV caller on all completed subgroups
```sh
cd /private/groups/patenlab/mira/centrolign/analysis/SVs_direct_pairwise

mkdir -p logs

sbatch /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/call_SVs_pairwise_all_chroms.sh /private/groups/patenlab/mira/centrolign/analysis/SVs_direct_pairwise/11172025_completed_subgroups.csv /private/groups/patenlab/mira/centrolign/analysis/SVs_direct_pairwise
```

Create csv file formatted as clade,chr,path_to_SV_beds
```sh
cd /private/groups/patenlab/mira/centrolign/analysis/SVs_direct_pairwise

echo "clade,chr,path_to_SV_beds" > 11172025_clade_chr_sv_beds.csv
cat /private/groups/patenlab/mira/centrolign/analysis/SVs_direct_pairwise/11172025_completed_subgroups.csv | while IFS=',' read -r clade cigars sample_lists ; do
  CHR=$(echo $clade | cut -f1 -d"_")
  echo $clade,$CHR,/private/groups/patenlab/mira/centrolign/analysis/SVs_direct_pairwise/${CHR}/SV_beds/${clade}/ >> 11172025_clade_chr_sv_beds.csv
done
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

#### Creating synteny plot with fedor's annotations overlaid

Convert fedor's SV beds to array coordinates
```sh
# get start coord of alpha sat sequence
grep chr12 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/HG00232.2_asat_arrays.bed
# HG00232#2#CM090029.1	34750702	37285438	chr12

grep chr12 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/HG01175.1_asat_arrays.bed
# HG01175#1#CM087931.1	34773369	37339781	chr12

# create HG00232.2 bed file
#  select Deletions
awk -F'\t' -v OFS='\t' -v n=34750702 '{ $2 -= n; $3 -= n }1' /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/fedor_horHap_SV_beds/HG00232.2_HG01175.1.bed | cut -f1-3,7,8 | awk -F'\t' -v OFS='\t' '{ print $0, "+", "0", "0", "0,0,0" }' | grep "D" > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_horHap_SVs.array_coords.bed

# create HG01175.1 bed file
#  select Insertions
awk -F'\t' -v OFS='\t' -v n=34773369 '{ $5 -= n; $6 -= n }1' /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/fedor_horHap_SV_beds/HG00232.2_HG01175.1.bed | cut -f4-6,7,8 | awk -F'\t' -v OFS='\t' '{ print $0, "+", "0", "0", "0,0,0" }' | grep "I" > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG01175.1_horHap_SVs.array_coords.bed
```

Run synteny plots
```sh
conda activate synteny

python /private/groups/migalab/juklucas/centrolign/chr12_test125/synteny_plot_bokeh.py   \
    --beds \
        /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_horHap_SVs.array_coords.bed \
        /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG01175.1_horHap_SVs.array_coords.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/induced_pairwise_cigars/pairwise_cigar_HG00232.2_HG01175.1.txt \
    --output /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/chr12_cenhap4_HG00232.2_HG01175.1_synteny.html \
    --show-mismatches \
    --web
```

Include horHap annotations as well.
```sh
awk -F'\t' -v OFS='\t' -v n=34750702 '{ $1="HG00232.2" ; $2 -= n; $3 -= n }1' /private/groups/migalab/fryabov/AS_annotation/cen12/horhap_annotation_align/bed/HG00232.2.bed | tail -n +4 > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_horhap_annotation_array_coords.bed

awk -F'\t' -v OFS='\t' -v n=34773369 '{ $1="HG01175.1" ; $2 -= n; $3 -= n }1' /private/groups/migalab/fryabov/AS_annotation/cen12/horhap_annotation_align/bed/HG01175.1.bed | tail -n +4 > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG01175.1_horhap_annotation_array_coords.bed

cat /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_horhap_annotation_array_coords.bed /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_horHap_SVs.array_coords.bed > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_horhap_SVs_annotations.bed

cat /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG01175.1_horhap_annotation_array_coords.bed /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG01175.1_horHap_SVs.array_coords.bed > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG01175.1_horhap_SVs_annotations.bed
```

Run synteny plots
```sh
conda activate synteny

python /private/groups/migalab/juklucas/centrolign/chr12_test125/synteny_plot_bokeh.py   \
    --beds \
        /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_horhap_SVs_annotations.bed \
        /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG01175.1_horhap_SVs_annotations.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/induced_pairwise_cigars/pairwise_cigar_HG00232.2_HG01175.1.txt \
    --output /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/chr12_cenhap4_HG00232.2_HG01175.1_synteny.with_horhaps.html \
    --show-mismatches \
    --web
```

### identify low divergence clades for chr 1 for Karen:

```sh
docker run -u `id -u`:`id -g` -v /private/groups:/private/groups/ \
    miramastoras/centromere_scripts:v0.1.4 \
    python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/identify_low_divergence_clades.py \
  --newick /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr1_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  --distance_csv /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/chr1_r2_QC_v2_centrolign_pairwise_distance.csv \
  --max_pairwise_dist 0.8 \
  --output_prefix /private/groups/patenlab/mira/chr1_low_divergence_clades
```

### Chromosome 12 case study

Stratifying SV calls by local identity

Make histogram of number of bases in parallelograms,traps,triangles overlapping with different local identity bins, at different pairwise distances

#### Intersect SV calls with local identity annotations

```sh
#!/bin/bash
#SBATCH --job-name=local_identity_intersect
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/call_SVs_%x.%j.log
#SBATCH --time=1:00:00
#SBATCH --array=[1]%1

# activate environment for bedtools
source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate base

SV_BED_PATHS=$1 # csv file where each line is a path to the SV bedPE file for each pair
OUTDIR=$2 # outdir for results
LOCAL_ID_IDX=$3 # local identity index file

SV_BED=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$SV_BED_PATHS")
REF_SMP=`echo $SV_BED | basename $SV_BED | cut -f1 -d"_"`
QUERY_SMP=`echo $SV_BED | basename $SV_BED | cut -f2 -d"_" | cut -f1-2 -d"."`

echo $SV_BED
echo $REF_SMP
echo $QUERY_SMP

# convert bedPE into bed file
LOCAL_FOLDER=/data/tmp/$(whoami)/${REF_SMP}_${QUERY_SMP}_tmp/
mkdir -p ${LOCAL_FOLDER}

cut -f1-3,7,8 $SV_BED > ${LOCAL_FOLDER}/all_SVs.bed
cut -f4-8 $SV_BED >> ${LOCAL_FOLDER}/all_SVs.bed

# sep by triangle, trap, parallelogram
awk '$5 < 0.1' ${LOCAL_FOLDER}/all_SVs.bed |  awk '$5 > 0' > ${LOCAL_FOLDER}/parallelograms.bed
awk '$5 > 0.1' ${LOCAL_FOLDER}/all_SVs.bed > ${LOCAL_FOLDER}/trapezoids.bed
awk '$5 == -1' ${LOCAL_FOLDER}/all_SVs.bed > ${LOCAL_FOLDER}/triangles.bed

# Get local identity bed file paths for each sample
REF_ID=`echo $REF_SMP | sed 's/\./#/g'`
QUERY_ID=`echo $QUERY_SMP | sed 's/\./#/g'`

REF_LOCAL_ID=`grep $REF_ID $LOCAL_ID_IDX | cut -f11 -d","`
QUERY_LOCAL_ID=`grep $REF_ID $LOCAL_ID_IDX | cut -f11 -d","`

# concatenate local identity tracks for both samples, convert contig name to match SV bed file
awk -v v="$REF_SMP" 'BEGIN { OFS="\t" } {$1=v; print}' ${REF_LOCAL_ID} > ${LOCAL_FOLDER}/local_identity_combined.bed
awk -v v="$QUERY_SMP" 'BEGIN { OFS="\t" } {$1=v; print}' ${QUERY_LOCAL_ID} >> ${LOCAL_FOLDER}/local_identity_combined.bed

# bedtools coverage to produce number of bases of each variant type overlapping each local identity window
bedtools coverage \
  -a ${LOCAL_FOLDER}/local_identity_combined.bed \
  -b ${LOCAL_FOLDER}/parallelograms.bed \
  > ${LOCAL_FOLDER}/local_iden_par.bed

bedtools coverage \
  -a ${LOCAL_FOLDER}/local_identity_combined.bed \
  -b ${LOCAL_FOLDER}/trapezoids.bed \
  > ${LOCAL_FOLDER}/local_iden_traps.bed

bedtools coverage \
  -a ${LOCAL_FOLDER}/local_identity_combined.bed \
  -b ${LOCAL_FOLDER}/triangles.bed \
  > ${LOCAL_FOLDER}/local_iden_tri.bed

cut -f1-9,11 ${LOCAL_FOLDER}/local_iden_par.bed > ${LOCAL_FOLDER}/par_counts.bed
cut -f11 ${LOCAL_FOLDER}/local_iden_traps.bed > ${LOCAL_FOLDER}/traps_counts.bed
cut -f11 ${LOCAL_FOLDER}/local_iden_tri.bed > ${LOCAL_FOLDER}/tri_counts.bed

paste -d"\t" ${LOCAL_FOLDER}/par_counts.bed ${LOCAL_FOLDER}/traps_counts.bed ${LOCAL_FOLDER}/tri_counts.bed > ${OUTDIR}/${REF_SMP}_${QUERY_SMP}.local_identity_SV_counts.bed

# clean up
rm -rf $LOCAL_FOLDER
```

```sh
# prepare input file listing out all bedPEs
cd /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/SV_beds/chr12_subgroup0/
ls | while read line ; do realpath $line ; done > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/SV_beds/local_identity/chr12_subgroup0_beds.txt

cd /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/SV_beds/chr12_subgroup1/
ls | while read line ; do realpath $line ; done > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/SV_beds/local_identity/chr12_subgroup1_beds.txt

# set up directory
mkdir -p /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/SV_beds/local_identity/chr12_subgroup_0
mkdir -p /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/SV_beds/local_identity/chr12_subgroup_1

cd /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/SV_beds/local_identity/
mkdir -p logs

# submit script to get count of bases within each local identity window for each sample pair
# output bed files will have # of bases from parallelograms, trapezoids, triangles as last 3 rows in that order
sbatch local_identity_intersect.sh \
    /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/SV_beds/local_identity/chr12_subgroup0_beds.txt \
    /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/SV_beds/local_identity/chr12_subgroup_0 \
    /private/groups/migalab/juklucas/censat_paper/local_identity/2025_11_13/chr12/chr12_local_id_launch_outputs.csv

sbatch local_identity_intersect.sh \
    /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/SV_beds/local_identity/chr12_subgroup1_beds.txt \
    /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/SV_beds/local_identity/chr12_subgroup_1 \
    /private/groups/migalab/juklucas/censat_paper/local_identity/2025_11_13/chr12/chr12_local_id_launch_outputs.csv
```

1. convert bedpe into bed file, then separate by triangle, trap, parallelogram.

2. for each pair, merge local identity tracks for both samples

3. bedtools intersect, return the identity track with the count of bases overlapping each window

4. Read in bed files into pandas df

5. Plot histogram of number of bases in parallelograms,traps,triangles overlapping with different local identity bins


Plotting SV length by local identity

Alignment breakpoints by HORHap for cenhap4
count of SV breakpoints by HORHAP broken by the categories
