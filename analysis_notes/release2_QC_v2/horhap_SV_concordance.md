## This notebook looks at concordance between HorHap derived SVs and centrolign SVs

### 1. Comparison using percent overlap

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

##### Creating synteny plot with fedor's annotations overlaid

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

Because INDEL representation is likely the cause of the low concordance, we need a slop added to our SV windows for concordance  

### 2. Measuring concordance using a ratio within windows

#### Taking example alignment between HG01175.1 and HG00232.2 to select window sizes

Expand HorHaps by 100, 1000 bp
```sh
source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate base

awk -v OFS="\t" '{print $1,$2-100,$3+100,$4,$5,$6,$7,$8,$9}' /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_horHap_SVs.array_coords.bed > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_horHap_SVs.array_coords.slop100.bed

awk -v OFS="\t" '{print $1,$2-100,$3+100,$4,$5,$6,$7,$8,$9}' /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG01175.1_horHap_SVs.array_coords.bed > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG01175.1_horHap_SVs.array_coords.slop100.bed

awk -v OFS="\t" '{print $1,$2-1000,$3+1000,$4,$5,$6,$7,$8,$9}' /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_horHap_SVs.array_coords.bed > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_horHap_SVs.array_coords.slop1000.bed

awk -v OFS="\t" '{print $1,$2-1000,$3+1000,$4,$5,$6,$7,$8,$9}' /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG01175.1_horHap_SVs.array_coords.bed > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG01175.1_horHap_SVs.array_coords.slop1000.bed
```

Take all centrolign SVs and expand coordinates.
```sh
grep "I" /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/SV_beds/chr12_subgroup1/HG00232.2_HG01175.1.bed | cut -f 4,5,6 | awk -v OFS="\t" '{print $1,$2-1000,$3+1000}' > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_HG01175.1.centrolign_SVs_slop1000.bed

grep "D" /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/SV_beds/chr12_subgroup1/HG00232.2_HG01175.1.bed | cut -f 1,2,3 | awk -v OFS="\t" '{print $1,$2-1000,$3+1000}' >> /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_HG01175.1.centrolign_SVs_slop1000.bed
```


Merge Horhap and centrolign SV windows
```sh
cat /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_horHap_SVs.array_coords.slop1000.bed /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG01175.1_horHap_SVs.array_coords.slop1000.bed | cut -f1-3 > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_HG01175.1.horhap_SVs_slop1000.bed

cat /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_HG01175.1.centrolign_SVs_slop1000.bed /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_HG01175.1.horhap_SVs_slop1000.bed | bedtools sort -i - | bedtools merge -i - > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_HG01175.1.1kb_windows.bed

grep HG00232.2 /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_HG01175.1.1kb_windows.bed | awk -F'\t' -v OFS='\t' '{ print $0, "SV","0","+", "0", "0", "0,0,0" }' > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2.1kb_windows.bed

grep HG01175.1 /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_HG01175.1.1kb_windows.bed | awk -F'\t' -v OFS='\t' '{ print $0, "SV","0","+", "0", "0", "0,0,0" }' > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG01175.1.1kb_windows.bed
```

```sh
grep "I" /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/SV_beds/chr12_subgroup1/HG00232.2_HG01175.1.bed | cut -f 4,5,6 | awk -v OFS="\t" '{print $1,$2-100,$3+100}' > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_HG01175.1.centrolign_SVs_slop100.bed

grep "D" /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/SV_beds/chr12_subgroup1/HG00232.2_HG01175.1.bed | cut -f 1,2,3 | awk -v OFS="\t" '{print $1,$2-100,$3+100}' >> /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_HG01175.1.centrolign_SVs_slop100.bed

cat /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_horHap_SVs.array_coords.slop100.bed /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG01175.1_horHap_SVs.array_coords.slop100.bed | cut -f1-3 > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_HG01175.1.horhap_SVs_slop100.bed

cat /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_HG01175.1.centrolign_SVs_slop100.bed /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_HG01175.1.horhap_SVs_slop100.bed | bedtools sort -i - | bedtools merge -i - > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_HG01175.1.100bp_windows.bed

grep HG00232.2 /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_HG01175.1.100bp_windows.bed | awk -F'\t' -v OFS='\t' '{ print $0, "SV","0","+", "0", "0", "0,0,0" }' > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2.100bp_windows.bed

grep HG01175.1 /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_HG01175.1.100bp_windows.bed | awk -F'\t' -v OFS='\t' '{ print $0, "SV","0","+", "0", "0", "0,0,0" }' > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG01175.1.100bp_windows.bed
```

Visualize overlap alongside new concordance windows
```sh
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/plot_bed_tracks.py \
  --region HG00232.2:0-2534736 \
  --beds /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2.d.bed /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_horHap_SVs.array_coords.bed /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2.100bp_windows.bed /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2.1kb_windows.bed \
  --labels centrolign horhap 100bp 1kb \
  --output /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00232.2_tracks.png
```

![plots](plots/HG00232.2_tracks.png)


#### Implementing windows of 100bp and 1000 bp for all cenhap 4 samples, calculating ratio of concordance


```sh
#!/bin/bash
#SBATCH --job-name=horhap_SV_concordance
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/call_SVs_%x.%j.log
#SBATCH --time=12:00:00

source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate base

OUTDIR=/private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/chr12/cenhap4/
mkdir -p $OUTDIR

HORHAP_BEDPE_DIR=/private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/fedor_horHap_SV_beds
CENTROLIGN_BEDS=/private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/SV_beds_asm_coords

cd ${CENTROLIGN_BEDS}

LOCAL_FOLDER=/data/tmp/$(whoami)/chr12_cenhap4/
mkdir -p ${LOCAL_FOLDER}

ls *.bed | while read -r bed ; do

  # for each bed get ref and query sample ID
  REF_SMP=`echo $bed | basename $bed | cut -f1 -d"_"`
  QUERY_SMP=`echo $bed | basename $bed | cut -f2 -d"_" | cut -f1-2 -d"."`

  # make sure horhap bed exists for the pair, and figure out which is ref and which is query
  REV="${HORHAP_BEDPE_DIR}/${QUERY_SMP}_${REF_SMP}.bed"
  MATCH="${HORHAP_BEDPE_DIR}/${bed}"

  if [[ -f "$MATCH" ]]; then
        HORHAP_BEDPE="$MATCH"
    elif [[ -f "$REV" ]]; then
        HORHAP_BEDPE="$REV"
    else
        echo "Skipping ${REF_SMP}, ${QUERY_SMP}, missing beds"
        continue
    fi

  # convert centrolign and horhap bed into bed pe files
  CEN_BEDPE=${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.triangles.bedpe
  awk '$8 == -1' $bed > $CEN_BEDPE

  CEN_BED=${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.triangles.bed
  grep "I" ${CEN_BEDPE} | cut -f 4,5,6 > $CEN_BED
  grep "D" ${CEN_BEDPE} | cut -f 1,2,3 >> $CEN_BED

  HORHAP_BED=${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.horhap.bed
  grep "I" ${HORHAP_BEDPE} | cut -f 4,5,6 > $HORHAP_BED
  grep "D" ${HORHAP_BEDPE} | cut -f 1,2,3 >> $HORHAP_BED

  ### Create 100bp windows around both sets of SVs ###
  # Separate centrolign SVs by insertions and deletions, add 100 bp slop
  grep "I" ${CEN_BEDPE} | cut -f 4,5,6 | awk -v OFS="\t" '{print $1,$2-100,$3+100}' > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.100bp.bed
  grep "D" ${CEN_BEDPE} | cut -f 1,2,3 | awk -v OFS="\t" '{print $1,$2-100,$3+100}' >> ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.100bp.bed

  grep "I" ${HORHAP_BEDPE} | cut -f 4,5,6 | awk -v OFS="\t" '{print $1,$2-100,$3+100}' >> ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.100bp.bed
  grep "D" ${HORHAP_BEDPE} | cut -f 1,2,3 | awk -v OFS="\t" '{print $1,$2-100,$3+100}' >> ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.100bp.bed

  # Merge windows for horhap and centrolign SVs
  bedtools sort -i ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.100bp.bed | bedtools merge -i - > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.100bp.windows.bed

  # Use bedtools coverage to calculate the number of bases in each window
  bedtools coverage \
    -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.100bp.windows.bed \
    -b ${CEN_BED} \
    > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.centrolign_cov.bed

  # Use bedtools coverage to calculate the number of bases in each window
  bedtools coverage \
    -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.100bp.windows.bed \
    -b ${HORHAP_BED} \
    > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.horhap_cov.bed

  # Combine coverage columns
  cut -f1-3,5 ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.centrolign_cov.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.centrolign_counts.bed

  cut -f5 ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.horhap_cov.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.horhap_counts.bed

  paste -d"\t" ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.centrolign_counts.bed ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.horhap_counts.bed > ${OUTDIR}/${REF_SMP}_${QUERY_SMP}.triangles.100bp.bed
done

rm -rf ${LOCAL_FOLDER}/
```

100 bp, all SVs
```sh
#!/bin/bash
#SBATCH --job-name=horhap_SV_concordance_all_svs_100bp
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/call_SVs_%x.%j.log
#SBATCH --time=12:00:00

source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate base

OUTDIR=/private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/chr12/cenhap4/
mkdir -p $OUTDIR

HORHAP_BEDPE_DIR=/private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/fedor_horHap_SV_beds
CENTROLIGN_BEDS=/private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/SV_beds_asm_coords

cd ${CENTROLIGN_BEDS}

LOCAL_FOLDER=/data/tmp/$(whoami)/chr12_cenhap4_all_SVs_100bp/
mkdir -p ${LOCAL_FOLDER}

ls *.bed | while read -r bed ; do

  # for each bed get ref and query sample ID
  REF_SMP=`echo $bed | basename $bed | cut -f1 -d"_"`
  QUERY_SMP=`echo $bed | basename $bed | cut -f2 -d"_" | cut -f1-2 -d"."`

  # make sure horhap bed exists for the pair, and figure out which is ref and which is query
  REV="${HORHAP_BEDPE_DIR}/${QUERY_SMP}_${REF_SMP}.bed"
  MATCH="${HORHAP_BEDPE_DIR}/${bed}"

  if [[ -f "$MATCH" ]]; then
        HORHAP_BEDPE="$MATCH"
    elif [[ -f "$REV" ]]; then
        HORHAP_BEDPE="$REV"
    else
        echo "Skipping ${REF_SMP}, ${QUERY_SMP}, missing beds"
        continue
    fi

  # use all SVs, not just triangles
  # convert cenbed to bed pe
  CEN_BEDPE=${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_SVs.bedpe
  cp $bed $CEN_BEDPE

  CEN_BED=${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_SVs.bed
  grep "I" ${CEN_BEDPE} | cut -f 4,5,6 > $CEN_BED
  grep "D" ${CEN_BEDPE} | cut -f 1,2,3 >> $CEN_BED

  HORHAP_BED=${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.horhap.bed
  grep "I" ${HORHAP_BEDPE} | cut -f 4,5,6 > $HORHAP_BED
  grep "D" ${HORHAP_BEDPE} | cut -f 1,2,3 >> $HORHAP_BED

  ### Create 100bp windows around both sets of SVs ###
  # Separate centrolign SVs by insertions and deletions, add 100 bp slop
  grep "I" ${CEN_BEDPE} | cut -f 4,5,6 | awk -v OFS="\t" '{print $1,$2-100,$3+100}' > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.100bp.bed
  grep "D" ${CEN_BEDPE} | cut -f 1,2,3 | awk -v OFS="\t" '{print $1,$2-100,$3+100}' >> ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.100bp.bed

  grep "I" ${HORHAP_BEDPE} | cut -f 4,5,6 | awk -v OFS="\t" '{print $1,$2-100,$3+100}' >> ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.100bp.bed
  grep "D" ${HORHAP_BEDPE} | cut -f 1,2,3 | awk -v OFS="\t" '{print $1,$2-100,$3+100}' >> ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.100bp.bed

  # Merge windows for horhap and centrolign SVs
  bedtools sort -i ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.100bp.bed | bedtools merge -i - > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.100bp.windows.bed

  # Use bedtools coverage to calculate the number of bases in each window
  bedtools coverage \
    -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.100bp.windows.bed \
    -b ${CEN_BED} \
    > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.centrolign_cov.bed

  # Use bedtools coverage to calculate the number of bases in each window
  bedtools coverage \
    -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.100bp.windows.bed \
    -b ${HORHAP_BED} \
    > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.horhap_cov.bed

  # Combine coverage columns
  cut -f1-3,5 ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.centrolign_cov.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.centrolign_counts.bed

  cut -f5 ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.horhap_cov.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.horhap_counts.bed

  paste -d"\t" ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.centrolign_counts.bed ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.horhap_counts.bed > ${OUTDIR}/${REF_SMP}_${QUERY_SMP}.all_SVs.100bp.bed
done

rm -rf ${LOCAL_FOLDER}/
```

1000 bp
```sh
#!/bin/bash
#SBATCH --job-name=horhap_SV_concordance_1kb
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/call_SVs_%x.%j.log
#SBATCH --time=12:00:00

source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate base

OUTDIR=/private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/chr12/cenhap4/
mkdir -p $OUTDIR

HORHAP_BEDPE_DIR=/private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/fedor_horHap_SV_beds
CENTROLIGN_BEDS=/private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/SV_beds_asm_coords

cd ${CENTROLIGN_BEDS}

LOCAL_FOLDER=/data/tmp/$(whoami)/chr12_cenhap4_1000bp/
mkdir -p ${LOCAL_FOLDER}

ls *.bed | while read -r bed ; do

  # for each bed get ref and query sample ID
  REF_SMP=`echo $bed | basename $bed | cut -f1 -d"_"`
  QUERY_SMP=`echo $bed | basename $bed | cut -f2 -d"_" | cut -f1-2 -d"."`

  # make sure horhap bed exists for the pair, and figure out which is ref and which is query
  REV="${HORHAP_BEDPE_DIR}/${QUERY_SMP}_${REF_SMP}.bed"
  MATCH="${HORHAP_BEDPE_DIR}/${bed}"

  if [[ -f "$MATCH" ]]; then
        HORHAP_BEDPE="$MATCH"
    elif [[ -f "$REV" ]]; then
        HORHAP_BEDPE="$REV"
    else
        echo "Skipping ${REF_SMP}, ${QUERY_SMP}, missing beds"
        continue
    fi

  # convert centrolign and horhap bed into bed pe files
  CEN_BEDPE=${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.triangles.bedpe
  awk '$8 == -1' $bed > $CEN_BEDPE

  CEN_BED=${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.triangles.bed
  grep "I" ${CEN_BEDPE} | cut -f 4,5,6 > $CEN_BED
  grep "D" ${CEN_BEDPE} | cut -f 1,2,3 >> $CEN_BED

  HORHAP_BED=${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.horhap.bed
  grep "I" ${HORHAP_BEDPE} | cut -f 4,5,6 > $HORHAP_BED
  grep "D" ${HORHAP_BEDPE} | cut -f 1,2,3 >> $HORHAP_BED

  ### Create 1000bp windows around both sets of SVs ###
  # Separate centrolign SVs by insertions and deletions, add 1000 bp slop
  grep "I" ${CEN_BEDPE} | cut -f 4,5,6 | awk -v OFS="\t" '{print $1,$2-1000,$3+1000}' > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.1000bp.bed
  grep "D" ${CEN_BEDPE} | cut -f 1,2,3 | awk -v OFS="\t" '{print $1,$2-1000,$3+1000}' >> ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.1000bp.bed

  grep "I" ${HORHAP_BEDPE} | cut -f 4,5,6 | awk -v OFS="\t" '{print $1,$2-1000,$3+1000}' >> ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.1000bp.bed
  grep "D" ${HORHAP_BEDPE} | cut -f 1,2,3 | awk -v OFS="\t" '{print $1,$2-1000,$3+1000}' >> ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.1000bp.bed

  # Merge windows for horhap and centrolign SVs
  bedtools sort -i ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.1000bp.bed | bedtools merge -i - > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.1000bp.windows.bed

  # Use bedtools coverage to calculate the number of bases in each window
  bedtools coverage \
    -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.1000bp.windows.bed \
    -b ${CEN_BED} \
    > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.centrolign_cov.bed

  # Use bedtools coverage to calculate the number of bases in each window
  bedtools coverage \
    -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.1000bp.windows.bed \
    -b ${HORHAP_BED} \
    > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.horhap_cov.bed

  # Combine coverage columns
  cut -f1-3,5 ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.centrolign_cov.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.centrolign_counts.bed

  cut -f5 ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.horhap_cov.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.horhap_counts.bed

  paste -d"\t" ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.centrolign_counts.bed ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.horhap_counts.bed > ${OUTDIR}/${REF_SMP}_${QUERY_SMP}.triangles.1000bp.bed
done

rm -rf ${LOCAL_FOLDER}/
```

1000 bp, all SVs
```sh
#!/bin/bash
#SBATCH --job-name=horhap_SV_concordance_1kb_all_SVs
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/call_SVs_%x.%j.log
#SBATCH --time=12:00:00

source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate base

OUTDIR=/private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/chr12/cenhap4/
mkdir -p $OUTDIR

HORHAP_BEDPE_DIR=/private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/fedor_horHap_SV_beds
CENTROLIGN_BEDS=/private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/SV_beds_asm_coords

cd ${CENTROLIGN_BEDS}

LOCAL_FOLDER=/data/tmp/$(whoami)/chr12_cenhap4_1000bp_all_SVs/
mkdir -p ${LOCAL_FOLDER}

ls *.bed | while read -r bed ; do

  # for each bed get ref and query sample ID
  REF_SMP=`echo $bed | basename $bed | cut -f1 -d"_"`
  QUERY_SMP=`echo $bed | basename $bed | cut -f2 -d"_" | cut -f1-2 -d"."`

  # make sure horhap bed exists for the pair, and figure out which is ref and which is query
  REV="${HORHAP_BEDPE_DIR}/${QUERY_SMP}_${REF_SMP}.bed"
  MATCH="${HORHAP_BEDPE_DIR}/${bed}"

  if [[ -f "$MATCH" ]]; then
        HORHAP_BEDPE="$MATCH"
    elif [[ -f "$REV" ]]; then
        HORHAP_BEDPE="$REV"
    else
        echo "Skipping ${REF_SMP}, ${QUERY_SMP}, missing beds"
        continue
    fi

  # use all SVs, not just triangles
  # convert cenbed to bed pe
  CEN_BEDPE=${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_SVs.bedpe
  cp $bed $CEN_BEDPE

  CEN_BED=${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.all_SVs.bed
  grep "I" ${CEN_BEDPE} | cut -f 4,5,6 > $CEN_BED
  grep "D" ${CEN_BEDPE} | cut -f 1,2,3 >> $CEN_BED

  HORHAP_BED=${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.horhap.bed
  grep "I" ${HORHAP_BEDPE} | cut -f 4,5,6 > $HORHAP_BED
  grep "D" ${HORHAP_BEDPE} | cut -f 1,2,3 >> $HORHAP_BED

  ### Create 1000bp windows around both sets of SVs ###
  # Separate centrolign SVs by insertions and deletions, add 1000 bp slop
  grep "I" ${CEN_BEDPE} | cut -f 4,5,6 | awk -v OFS="\t" '{print $1,$2-1000,$3+1000}' > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.1000bp.bed
  grep "D" ${CEN_BEDPE} | cut -f 1,2,3 | awk -v OFS="\t" '{print $1,$2-1000,$3+1000}' >> ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.1000bp.bed

  grep "I" ${HORHAP_BEDPE} | cut -f 4,5,6 | awk -v OFS="\t" '{print $1,$2-1000,$3+1000}' >> ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.1000bp.bed
  grep "D" ${HORHAP_BEDPE} | cut -f 1,2,3 | awk -v OFS="\t" '{print $1,$2-1000,$3+1000}' >> ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.1000bp.bed

  # Merge windows for horhap and centrolign SVs
  bedtools sort -i ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.1000bp.bed | bedtools merge -i - > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.1000bp.windows.bed

  # Use bedtools coverage to calculate the number of bases in each window
  bedtools coverage \
    -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.1000bp.windows.bed \
    -b ${CEN_BED} \
    > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.centrolign_cov.bed

  # Use bedtools coverage to calculate the number of bases in each window
  bedtools coverage \
    -a ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.1000bp.windows.bed \
    -b ${HORHAP_BED} \
    > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.horhap_cov.bed

  # Combine coverage columns
  cut -f1-3,5 ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.centrolign_cov.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.centrolign_counts.bed

  cut -f5 ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.horhap_cov.bed > ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.horhap_counts.bed

  paste -d"\t" ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.centrolign_counts.bed ${LOCAL_FOLDER}/${REF_SMP}_${QUERY_SMP}.horhap_counts.bed > ${OUTDIR}/${REF_SMP}_${QUERY_SMP}.all_SVs.1000bp.bed
done

rm -rf ${LOCAL_FOLDER}/
```
### Plot synteny plots where concordance is the worst

Convert fedor's SV beds to array coordinates
```sh
# get start coord of alpha sat sequence
grep chr12 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/HG00290.2_asat_arrays.bed
# HG00290#2#CM090009.1	34720173	37437832	chr12

grep chr12 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/HG03710.1_asat_arrays.bed
# HG03710#1#CM086873.1	34650988	37479323	chr12

# create HG00290.2 bed file
#  select Deletions
awk -F'\t' -v OFS='\t' -v n=34720173 '{ $2 -= n; $3 -= n }1' /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/fedor_horHap_SV_beds/HG00290.2_HG03710.1.bed | cut -f1-3,7,8 | awk -F'\t' -v OFS='\t' '{ print $0, "+", "0", "0", "0,0,0" }' | grep "D" > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00290.2_horHap_SVs.array_coords.bed

# create HG03710.1 bed file
#  select Insertions
awk -F'\t' -v OFS='\t' -v n=34650988 '{ $5 -= n; $6 -= n }1' /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/fedor_horHap_SV_beds/HG00290.2_HG03710.1.bed | cut -f4-6,7,8 | awk -F'\t' -v OFS='\t' '{ print $0, "+", "0", "0", "0,0,0" }' | grep "I" > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG03710.1_horHap_SVs.array_coords.bed
```

Run synteny plots
```sh
conda activate synteny

python /private/groups/migalab/juklucas/centrolign/chr12_test125/synteny_plot_bokeh.py   \
    --beds \
        /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00290.2_horHap_SVs.array_coords.bed \
        /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG03710.1_horHap_SVs.array_coords.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/induced_pairwise_cigars/pairwise_cigar_HG00290.2_HG03710.1.txt \
    --output /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/chr12_cenhap4_HG00290.2_HG03710.1_synteny.html \
    --show-mismatches \
    --web
```

### Plot synteny plots where distance is 0.5 and concordance is great

Convert fedor's SV beds to array coordinates
```sh
# get start coord of alpha sat sequence
grep chr12 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/HG00329.2_asat_arrays.bed
# HG00329#2#CM094244.1	34789273	37579344	chr12

grep chr12 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/HG03654.2_asat_arrays.bed
# HG03654#2#CM086861.1	34643372	36478006	chr12

# create HG00329.2 bed file
#  select Deletions
awk -F'\t' -v OFS='\t' -v n=34789273 '{ $2 -= n; $3 -= n }1' /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/fedor_horHap_SV_beds/HG00329.2_HG03654.2.bed | cut -f1-3,7,8 | awk -F'\t' -v OFS='\t' '{ print $0, "+", "0", "0", "0,0,0" }' | grep "D" > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00329.2_horHap_SVs.array_coords.bed

# create HG03654.2 bed file
#  select Insertions
awk -F'\t' -v OFS='\t' -v n=34643372 '{ $5 -= n; $6 -= n }1' /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/fedor_horHap_SV_beds/HG00329.2_HG03654.2.bed | cut -f4-6,7,8 | awk -F'\t' -v OFS='\t' '{ print $0, "+", "0", "0", "0,0,0" }' | grep "I" > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG03654.2_horHap_SVs.array_coords.bed
```

Run synteny plots
```sh
conda activate synteny

python /private/groups/migalab/juklucas/centrolign/chr12_test125/synteny_plot_bokeh.py   \
    --beds \
        /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG00329.2_horHap_SVs.array_coords.bed \
        /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG03654.2_horHap_SVs.array_coords.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/induced_pairwise_cigars/pairwise_cigar_HG00329.2_HG03654.2.txt \
    --output /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/chr12_cenhap4_HG00329.2_HG03654.2_synteny.html \
    --show-mismatches \
    --web
```

### Plot synteny plots where distance is < 0.2 and concordance is great

Convert fedor's SV beds to array coordinates
```sh
# get start coord of alpha sat sequence
grep chr12 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/HG02698.1_asat_arrays.bed
# HG02698#1#CM086831.1	34776248	37256955	chr12

grep chr12 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/HG03784.2_asat_arrays.bed
# HG03784#2#CM094206.1	34781322	37439318	chr12

# create HG00329.2 bed file
#  select Deletions
awk -F'\t' -v OFS='\t' -v n=34776248 '{ $2 -= n; $3 -= n }1' /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/fedor_horHap_SV_beds/HG02698.1_HG03784.2.bed  | cut -f1-3,7,8 | awk -F'\t' -v OFS='\t' '{ print $0, "+", "0", "0", "0,0,0" }' | grep "D" > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG02698.1_horHap_SVs.array_coords.bed

# create HG03654.2 bed file
#  select Insertions
awk -F'\t' -v OFS='\t' -v n=34781322 '{ $5 -= n; $6 -= n }1' /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/fedor_horHap_SV_beds/HG02698.1_HG03784.2.bed | cut -f4-6,7,8 | awk -F'\t' -v OFS='\t' '{ print $0, "+", "0", "0", "0,0,0" }' | grep "I" > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG03784.2_horHap_SVs.array_coords.bed
```

Run synteny plots
```sh
conda activate synteny

python /private/groups/migalab/juklucas/centrolign/chr12_test125/synteny_plot_bokeh.py   \
    --beds \
        /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG02698.1_horHap_SVs.array_coords.bed \
        /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/HG03784.2_horHap_SVs.array_coords.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/induced_pairwise_cigars/pairwise_cigar_HG02698.1_HG03784.2.txt \
    --output /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/synteny_plots/chr12_cenhap4_HG02698.1_HG03784.2_synteny.html \
    --show-mismatches \
    --web
```

### New approach:

- For every centrolign SV, check if horhap SV exists within X bp up and downstream, that is within Y% of the size. If so, its correct, if not, incorrect
- repeat for horhap SV
- report the harmonic mean


```sh
# select only insertions from query, deletions from ref coords
grep "I" /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/SV_beds_asm_coords/HG02698.1_HG03784.2.bed | cut -f4,5,6,8 > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/test_f1/HG02698.1_HG03784.2.cen.bed

grep "D" /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/SV_beds_asm_coords/HG02698.1_HG03784.2.bed | cut -f1-3,8 >> /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/test_f1/HG02698.1_HG03784.2.cen.bed

awk '$4 == -1' /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/test_f1/HG02698.1_HG03784.2.cen.bed > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/test_f1/HG02698.1_HG03784.2.cen.tri.bed

grep "I" /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/fedor_horHap_SV_beds/HG02698.1_HG03784.2.bed | cut -f4,5,6,8 > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/test_f1/HG02698.1_HG03784.2.hor.bed

grep "D" /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/fedor_horHap_SV_beds/HG02698.1_HG03784.2.bed | cut -f1-3,8 >> /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/test_f1/HG02698.1_HG03784.2.hor.bed

python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/sv_compare.py \
    /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/test_f1/HG02698.1_HG03784.2.cen.bed \
    /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/test_f1/HG02698.1_HG03784.2.hor.bed \
    --max_dist 3000 --min_size_frac 0.5 --out_prefix /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/test_f1/results --bed9 --debug

#!/bin/sh
awk -F'\t' -v OFS='\t' -v n=34776248 '{ $2 -= n; $3 -= n }1' /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/test_f1/results_BED2.bed9 | grep HG02698.1 > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/test_f1/results_BED2_HG02698.1.bed9.array_coords.bed

awk -F'\t' -v OFS='\t' -v n=34781322 '{ $2 -= n; $3 -= n }1' /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/test_f1/results_BED2.bed9 | grep HG03784.2 > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/test_f1/results_BED2_HG03784.2.bed9.array_coords.bed
```

Synteny with matches/mismatches marked:
```sh
python /private/groups/migalab/juklucas/centrolign/chr12_test125/synteny_plot_bokeh.py   \
    --beds \
        /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/test_f1/results_BED2_HG02698.1.bed9.array_coords.bed \
        /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/test_f1/results_BED2_HG03784.2.bed9.array_coords.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/induced_pairwise_cigars/pairwise_cigar_HG02698.1_HG03784.2.txt \
    --output /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/test_f1/chr12_cenhap4_HG02698.1_HG03784.2_synteny.bed9.html \
    --show-mismatches \
    --web
```

Centrolign :

awk -F'\t' -v OFS='\t' -v n=34776248 '{ $2 -= n; $3 -= n }1' /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/test_f1/results_BED1.bed9 | grep HG02698.1 > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/test_f1/results_BED2_HG02698.1.bed9.array_coords.cen.bed

awk -F'\t' -v OFS='\t' -v n=34781322 '{ $2 -= n; $3 -= n }1' /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/test_f1/results_BED1.bed9 | grep HG03784.2 > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/test_f1/results_BED2_HG03784.2.bed9.array_coords.cen.bed
```

Synteny with matches/mismatches marked:
```sh
python /private/groups/migalab/juklucas/centrolign/chr12_test125/synteny_plot_bokeh.py   \
    --beds \
        /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/test_f1/results_BED2_HG02698.1.bed9.array_coords.cen.bed \
        /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/test_f1/results_BED2_HG03784.2.bed9.array_coords.cen.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/induced_pairwise_cigars/pairwise_cigar_HG02698.1_HG03784.2.txt \
    --output /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/test_f1/chr12_cenhap4_HG02698.1_HG03784.2_synteny.bed9.cen.html \
    --show-mismatches \
    --web
```


Run script for all cenhap4 samples


reverse all cigar strings in cenhap4 directory:
```sh

ls /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/induced_pairwise_cigars | while read line ; do
  python3 /private/groups/patenlab/mira/centrolign/github/censat_paper/scripts/centrolign_result_parsing/reverse_cigar.py \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/induced_pairwise_cigars/${line}
done
```

python3 /private/groups/patenlab/mira/centrolign/github/censat_paper/scripts/centrolign_result_parsing/reverse_cigar.py /private/groups/patenlab/mira/rev_cigar/

```sh
#!/bin/bash
#SBATCH --job-name=sv_compare
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --output=logs/sv_compare_%x.%j.log
#SBATCH --array=[1-1485]%100

BED_CSV=/private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/matched_beds.filt.csv

CENTROLIGN_BED=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$BED_CSV" | cut -f1 -d",")
HORHAP_BED=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$BED_CSV" | cut -f2 -d",")

SMP_PAIR=$(basename -s .bed $CENTROLIGN_BED )
SMP1="${SMP_PAIR%%_*}"
SMP2="${SMP_PAIR##*_}"  

echo $SMP_PAIR
echo $CENTROLIGN_BED
echo $HORHAP_BED

LOCAL_FOLDER=/data/tmp/$(whoami)/chr12_cenhap4_${SMP_PAIR}/
mkdir -p ${LOCAL_FOLDER}

OUTDIR=/private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/sv_compare/md3k_msf0.5
mkdir -p $OUTDIR/

# select only insertions from query, deletions from ref coords
grep "I" ${CENTROLIGN_BED} | cut -f4,5,6,8 > ${LOCAL_FOLDER}/${SMP_PAIR}.cen.bed

grep "D" ${CENTROLIGN_BED} | cut -f1-3,8 >> ${LOCAL_FOLDER}/${SMP_PAIR}.cen.bed

grep "I" ${HORHAP_BED} | cut -f4,5,6,8 > ${LOCAL_FOLDER}/${SMP_PAIR}.hor.bed

grep "D" ${HORHAP_BED} | cut -f1-3,8 >> ${LOCAL_FOLDER}/${SMP_PAIR}.hor.bed

# Run SV comparison
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/sv_compare.py \
    ${LOCAL_FOLDER}/${SMP_PAIR}.cen.bed \
    ${LOCAL_FOLDER}/${SMP_PAIR}.hor.bed \
    --max_dist 3000 --min_size_frac 0.5 \
    --out_prefix ${OUTDIR}/${SMP_PAIR} --bed9 --debug

# Auto-generate synteny plots showing cenhap and horhap SV qualities
# get start coord of alpha sat sequence
SMP1_START=`grep chr12 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/${SMP1}_asat_arrays.bed | cut -f2`

SMP2_START=`grep chr12 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/${SMP2}_asat_arrays.bed | cut -f2`

echo $SMP1_START
echo $SMP2_START

# HORHAP
awk -F'\t' -v OFS='\t' -v n=${SMP1_START} '{ $2 -= n; $3 -= n }1' ${OUTDIR}/${SMP_PAIR}_BED2.bed9 | grep ${SMP1} > ${LOCAL_FOLDER}/${SMP1}_sv_horhap.bed

awk -F'\t' -v OFS='\t' -v n=${SMP2_START} '{ $2 -= n; $3 -= n }1' ${OUTDIR}/${SMP_PAIR}_BED2.bed9 | grep ${SMP2} > ${LOCAL_FOLDER}/${SMP2}_sv_horhap.bed

# CENTROLIGN
awk -F'\t' -v OFS='\t' -v n=${SMP1_START} '{ $2 -= n; $3 -= n }1' ${OUTDIR}/${SMP_PAIR}_BED1.bed9 | grep ${SMP1} > ${LOCAL_FOLDER}/${SMP1}_sv_cen.bed

awk -F'\t' -v OFS='\t' -v n=${SMP2_START} '{ $2 -= n; $3 -= n }1' ${OUTDIR}/${SMP_PAIR}_BED1.bed9 | grep ${SMP2} > ${LOCAL_FOLDER}/${SMP2}_sv_cen.bed

# Generate synteny plots
source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate synteny

mkdir -p ${OUTDIR}/synteny_plots
# horhap view
python /private/groups/migalab/juklucas/centrolign/chr12_test125/synteny_plot_bokeh.py   \
    --beds \
        ${LOCAL_FOLDER}/${SMP1}_sv_horhap.bed \
        ${LOCAL_FOLDER}/${SMP2}_sv_horhap.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/induced_pairwise_cigars/pairwise_cigar_${SMP1}_${SMP2}.txt \
    --output ${OUTDIR}/synteny_plots/md3k_msf0.5_${SMP1}_${SMP2}_synteny.horhap.html \
    --show-mismatches \
    --web

# centrolign view
python /private/groups/migalab/juklucas/centrolign/chr12_test125/synteny_plot_bokeh.py   \
    --beds \
        ${LOCAL_FOLDER}/${SMP1}_sv_cen.bed \
        ${LOCAL_FOLDER}/${SMP2}_sv_cen.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/induced_pairwise_cigars/pairwise_cigar_${SMP1}_${SMP2}.txt \
    --output ${OUTDIR}/synteny_plots/md3k_msf0.5_${SMP1}_${SMP2}_synteny.cen.html \
    --show-mismatches \
    --web

rm -rf ${LOCAL_FOLDER}
```

Submitting script
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/sv_compare
cd /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/sv_compare

mkdir md3k_msf0.5
mkdir -p logs

sbatch sv_compare.sh
```

Look at synteny plots for samples with the worst F1 in bins 0-0.2

Example 1
```sh
# get start coord of alpha sat sequence
grep chr12 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/HG00290.2_asat_arrays.bed
# HG00290#2#CM090009.1	34720173	37437832	chr12

grep chr12 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/HG01940.1_asat_arrays.bed
# HG01940#1#CM088637.1	34729153	37567028	chr12

awk -F'\t' -v OFS='\t' -v n=34720173 '{ $2 -= n; $3 -= n }1' /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/sv_compare/md3k_msf0.5/HG00290.2_HG01940.1_BED2.bed9 | grep HG00290.2 > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/sv_compare/md3k_msf0.5/HG00290.2_sv_HG01940.1_horhap.bed

awk -F'\t' -v OFS='\t' -v n=34729153 '{ $2 -= n; $3 -= n }1' /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/sv_compare/md3k_msf0.5/HG00290.2_HG01940.1_BED2.bed9 | grep HG01940.1 > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/sv_compare/md3k_msf0.5/HG01940.1_sv_HG00290.2_horhap.bed

python /private/groups/migalab/juklucas/centrolign/chr12_test125/synteny_plot_bokeh.py   \
    --beds \
        /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/sv_compare/md3k_msf0.5/HG00290.2_sv_HG01940.1_horhap.bed \
         /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/sv_compare/md3k_msf0.5/HG01940.1_sv_HG00290.2_horhap.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/induced_pairwise_cigars/pairwise_cigar_HG00290.2_HG01940.1.txt \
    --output /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/sv_compare/synteny_plots/md3k_msf0.5_HG00290.2_HG01940.1_synteny.horhap.html \
    --show-mismatches \
    --web
```
Example 2
```sh
# get start coord of alpha sat sequence
grep chr12 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/HG00290.2_asat_arrays.bed
# HG00290#2#CM090009.1	34720173	37437832	chr12

grep chr12 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/HG01940.1_asat_arrays.bed
# HG01940#1#CM088637.1	34729153	37567028	chr12

awk -F'\t' -v OFS='\t' -v n=34720173 '{ $2 -= n; $3 -= n }1' /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/sv_compare/md3k_msf0.5/HG00290.2_HG01940.1_BED2.bed9 | grep HG00290.2 > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/sv_compare/md3k_msf0.5/HG00290.2_sv_HG01940.1_horhap.bed

awk -F'\t' -v OFS='\t' -v n=34729153 '{ $2 -= n; $3 -= n }1' /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/sv_compare/md3k_msf0.5/HG00290.2_HG01940.1_BED2.bed9 | grep HG01940.1 > /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/sv_compare/md3k_msf0.5/HG01940.1_sv_HG00290.2_horhap.bed

python /private/groups/migalab/juklucas/centrolign/chr12_test125/synteny_plot_bokeh.py   \
    --beds \
        /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/sv_compare/md3k_msf0.5/HG00290.2_sv_HG01940.1_horhap.bed \
         /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/sv_compare/md3k_msf0.5/HG01940.1_sv_HG00290.2_horhap.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/induced_pairwise_cigars/pairwise_cigar_HG00290.2_HG01940.1.txt \
    --output /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/sv_compare/synteny_plots/md3k_msf0.5_HG00290.2_HG01940.1_synteny.horhap.html \
    --show-mismatches \
    --web
```


### Get list of all sample pairs with distances < 0.4 for all chroms

```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices

ls | while read line ; do
    awk -F',' '$3 <= 0.4' $line > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices_lt0.4/$line
done

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices

ls | while read line ; do
    awk -F',' '$3 <= 0.2' $line > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices_lt0.2/$line
done
```
### Chr X horHap SV concordance


```sh
#!/bin/bash
#SBATCH --job-name=sv_compare_chrX
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --output=logs/sv_compare_%x.%j.log
#SBATCH --array=[1-6991]%100

BED_CSV=/private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/f1_scores/chrX/matched_beds.filt.csv

CENTROLIGN_BED=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$BED_CSV" | cut -f1 -d",")
HORHAP_BED=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$BED_CSV" | cut -f2 -d",")

SMP_PAIR=$(basename -s .bed $CENTROLIGN_BED )
SMP1="${SMP_PAIR%%_*}"
SMP2="${SMP_PAIR##*_}"  

echo $SMP_PAIR
echo $CENTROLIGN_BED
echo $HORHAP_BED

LOCAL_FOLDER=/data/tmp/$(whoami)/chrX_${SMP_PAIR}/
mkdir -p ${LOCAL_FOLDER}

OUTDIR=/private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/f1_scores/chrX/sv_compare/md3k_msf0.5
mkdir -p $OUTDIR/

# select only insertions from query, deletions from ref coords
grep "I" ${CENTROLIGN_BED} | cut -f4,5,6,8 > ${LOCAL_FOLDER}/${SMP_PAIR}.cen.bed

grep "D" ${CENTROLIGN_BED} | cut -f1-3,8 >> ${LOCAL_FOLDER}/${SMP_PAIR}.cen.bed

grep "I" ${HORHAP_BED} | cut -f4,5,6,8 > ${LOCAL_FOLDER}/${SMP_PAIR}.hor.bed

grep "D" ${HORHAP_BED} | cut -f1-3,8 >> ${LOCAL_FOLDER}/${SMP_PAIR}.hor.bed

# Run SV comparison
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/sv_compare.py \
    ${LOCAL_FOLDER}/${SMP_PAIR}.cen.bed \
    ${LOCAL_FOLDER}/${SMP_PAIR}.hor.bed \
    --max_dist 3000 --min_size_frac 0.5 \
    --out_prefix ${OUTDIR}/${SMP_PAIR} --bed9 --debug

# Auto-generate synteny plots showing cenhap and horhap SV qualities
# get start coord of alpha sat sequence
SMP1_START=`grep chrX /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/${SMP1}_asat_arrays.bed | cut -f2`

SMP2_START=`grep chrX /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/${SMP2}_asat_arrays.bed | cut -f2`

echo $SMP1_START
echo $SMP2_START

# HORHAP
awk -F'\t' -v OFS='\t' -v n=${SMP1_START} '{ $2 -= n; $3 -= n }1' ${OUTDIR}/${SMP_PAIR}_BED2.bed9 | grep ${SMP1} > ${LOCAL_FOLDER}/${SMP1}_sv_horhap.bed

awk -F'\t' -v OFS='\t' -v n=${SMP2_START} '{ $2 -= n; $3 -= n }1' ${OUTDIR}/${SMP_PAIR}_BED2.bed9 | grep ${SMP2} > ${LOCAL_FOLDER}/${SMP2}_sv_horhap.bed

# CENTROLIGN
awk -F'\t' -v OFS='\t' -v n=${SMP1_START} '{ $2 -= n; $3 -= n }1' ${OUTDIR}/${SMP_PAIR}_BED1.bed9 | grep ${SMP1} > ${LOCAL_FOLDER}/${SMP1}_sv_cen.bed

awk -F'\t' -v OFS='\t' -v n=${SMP2_START} '{ $2 -= n; $3 -= n }1' ${OUTDIR}/${SMP_PAIR}_BED1.bed9 | grep ${SMP2} > ${LOCAL_FOLDER}/${SMP2}_sv_cen.bed

# Generate synteny plots
source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate synteny

mkdir -p ${OUTDIR}/synteny_plots
# horhap view
python /private/groups/migalab/juklucas/centrolign/chr12_test125/synteny_plot_bokeh.py   \
    --beds \
        ${LOCAL_FOLDER}/${SMP1}_sv_horhap.bed \
        ${LOCAL_FOLDER}/${SMP2}_sv_horhap.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrX/induced_pairwise_cigars/pairwise_cigar_${SMP1}_${SMP2}.txt \
    --output ${OUTDIR}/synteny_plots/md3k_msf0.5_${SMP1}_${SMP2}_synteny.horhap.html \
    --show-mismatches \
    --web

# centrolign view
python /private/groups/migalab/juklucas/centrolign/chr12_test125/synteny_plot_bokeh.py   \
    --beds \
        ${LOCAL_FOLDER}/${SMP1}_sv_cen.bed \
        ${LOCAL_FOLDER}/${SMP2}_sv_cen.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrX/induced_pairwise_cigars/pairwise_cigar_${SMP1}_${SMP2}.txt \
    --output ${OUTDIR}/synteny_plots/md3k_msf0.5_${SMP1}_${SMP2}_synteny.cen.html \
    --show-mismatches \
    --web

rm -rf ${LOCAL_FOLDER}
```

Submitting script
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/sv_compare
cd /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/sv_compare

mkdir md3k_msf0.5
mkdir -p logs

sbatch sv_compare.sh
```

#### SV compare 2 - windows

Submitting script
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/sv_compare2
cd /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/sv_compare2

mkdir md3k_msf0.5
mkdir -p logs

sbatch sv_compare.sh
```

```sh
#!/bin/bash
#SBATCH --job-name=sv_compare
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --output=logs/sv_compare_%x.%j.log
#SBATCH --array=[1-1485]%100

BED_CSV=/private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr12/cenHap4_benchmarking_HorHaps/matched_beds.filt.csv

CENTROLIGN_BED=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$BED_CSV" | cut -f1 -d",")
HORHAP_BED=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$BED_CSV" | cut -f2 -d",")

SMP_PAIR=$(basename -s .bed $CENTROLIGN_BED )
SMP1="${SMP_PAIR%%_*}"
SMP2="${SMP_PAIR##*_}"  

echo $SMP_PAIR
echo $CENTROLIGN_BED
echo $HORHAP_BED

LOCAL_FOLDER=/data/tmp/$(whoami)/chr12_cenhap4_${SMP_PAIR}/
mkdir -p ${LOCAL_FOLDER}

OUTDIR=/private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/f1_scores/chr12/sv_compare2/md3k_msf0.5
mkdir -p $OUTDIR/

# select only insertions from query, deletions from ref coords
grep "I" ${CENTROLIGN_BED} | cut -f4,5,6,8 > ${LOCAL_FOLDER}/${SMP_PAIR}.cen.bed

grep "D" ${CENTROLIGN_BED} | cut -f1-3,8 >> ${LOCAL_FOLDER}/${SMP_PAIR}.cen.bed

grep "I" ${HORHAP_BED} | cut -f4,5,6,8 > ${LOCAL_FOLDER}/${SMP_PAIR}.hor.bed

grep "D" ${HORHAP_BED} | cut -f1-3,8 >> ${LOCAL_FOLDER}/${SMP_PAIR}.hor.bed

# Run SV comparison
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/sv_compare2.py \
    ${LOCAL_FOLDER}/${SMP_PAIR}.cen.bed \
    ${LOCAL_FOLDER}/${SMP_PAIR}.hor.bed \
    --max_dist 3000 --min_size_frac 0.5 \
    --out_prefix ${OUTDIR}/${SMP_PAIR} --bed9 --debug

# Auto-generate synteny plots showing cenhap and horhap SV qualities
# get start coord of alpha sat sequence
SMP1_START=`grep chr12 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/${SMP1}_asat_arrays.bed | cut -f2`

SMP2_START=`grep chr12 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/${SMP2}_asat_arrays.bed | cut -f2`

echo $SMP1_START
echo $SMP2_START

# HORHAP
awk -F'\t' -v OFS='\t' -v n=${SMP1_START} '{ $2 -= n; $3 -= n }1' ${OUTDIR}/${SMP_PAIR}_BED2.bed9 | grep ${SMP1} > ${LOCAL_FOLDER}/${SMP1}_sv_horhap.bed

awk -F'\t' -v OFS='\t' -v n=${SMP2_START} '{ $2 -= n; $3 -= n }1' ${OUTDIR}/${SMP_PAIR}_BED2.bed9 | grep ${SMP2} > ${LOCAL_FOLDER}/${SMP2}_sv_horhap.bed

# CENTROLIGN
awk -F'\t' -v OFS='\t' -v n=${SMP1_START} '{ $2 -= n; $3 -= n }1' ${OUTDIR}/${SMP_PAIR}_BED1.bed9 | grep ${SMP1} > ${LOCAL_FOLDER}/${SMP1}_sv_cen.bed

awk -F'\t' -v OFS='\t' -v n=${SMP2_START} '{ $2 -= n; $3 -= n }1' ${OUTDIR}/${SMP_PAIR}_BED1.bed9 | grep ${SMP2} > ${LOCAL_FOLDER}/${SMP2}_sv_cen.bed

# Generate synteny plots
source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate synteny

mkdir -p ${OUTDIR}/synteny_plots
# horhap view
python /private/groups/migalab/juklucas/centrolign/chr12_test125/synteny_plot_bokeh.py   \
    --beds \
        ${LOCAL_FOLDER}/${SMP1}_sv_horhap.bed \
        ${LOCAL_FOLDER}/${SMP2}_sv_horhap.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/induced_pairwise_cigars/pairwise_cigar_${SMP1}_${SMP2}.txt \
    --output ${OUTDIR}/synteny_plots/md3k_msf0.5_${SMP1}_${SMP2}_synteny.horhap.html \
    --show-mismatches \
    --web

# centrolign view
python /private/groups/migalab/juklucas/centrolign/chr12_test125/synteny_plot_bokeh.py   \
    --beds \
        ${LOCAL_FOLDER}/${SMP1}_sv_cen.bed \
        ${LOCAL_FOLDER}/${SMP2}_sv_cen.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/induced_pairwise_cigars/pairwise_cigar_${SMP1}_${SMP2}.txt \
    --output ${OUTDIR}/synteny_plots/md3k_msf0.5_${SMP1}_${SMP2}_synteny.cen.html \
    --show-mismatches \
    --web

rm -rf ${LOCAL_FOLDER}
```


Chr X SVcompare2

```sh
#!/bin/bash
#SBATCH --job-name=sv_compare2_chrX
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --output=logs/sv_compare_%x.%j.log
#SBATCH --array=[1-6991]%100

BED_CSV=/private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/f1_scores/chrX/matched_beds.filt.csv

CENTROLIGN_BED=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$BED_CSV" | cut -f1 -d",")
HORHAP_BED=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$BED_CSV" | cut -f2 -d",")

SMP_PAIR=$(basename -s .bed $CENTROLIGN_BED )
SMP1="${SMP_PAIR%%_*}"
SMP2="${SMP_PAIR##*_}"  

echo $SMP_PAIR
echo $CENTROLIGN_BED
echo $HORHAP_BED

LOCAL_FOLDER=/data/tmp/$(whoami)/chrX_${SMP_PAIR}/
mkdir -p ${LOCAL_FOLDER}

OUTDIR=/private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/f1_scores/chrX/sv_compare2/md3k_msf0.5
mkdir -p $OUTDIR/

# select only insertions from query, deletions from ref coords
grep "I" ${CENTROLIGN_BED} | cut -f4,5,6,8 > ${LOCAL_FOLDER}/${SMP_PAIR}.cen.bed

grep "D" ${CENTROLIGN_BED} | cut -f1-3,8 >> ${LOCAL_FOLDER}/${SMP_PAIR}.cen.bed

grep "I" ${HORHAP_BED} | cut -f4,5,6,8 > ${LOCAL_FOLDER}/${SMP_PAIR}.hor.bed

grep "D" ${HORHAP_BED} | cut -f1-3,8 >> ${LOCAL_FOLDER}/${SMP_PAIR}.hor.bed

# Run SV comparison
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/sv_comparev2.py \
    ${LOCAL_FOLDER}/${SMP_PAIR}.cen.bed \
    ${LOCAL_FOLDER}/${SMP_PAIR}.hor.bed \
    --max_dist 3000 --min_size_frac 0.5 \
    --out_prefix ${OUTDIR}/${SMP_PAIR} --bed9 --debug

# Auto-generate synteny plots showing cenhap and horhap SV qualities
# get start coord of alpha sat sequence
SMP1_START=`grep chrX /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/${SMP1}_asat_arrays.bed | cut -f2`

SMP2_START=`grep chrX /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/${SMP2}_asat_arrays.bed | cut -f2`

echo $SMP1_START
echo $SMP2_START

# HORHAP
awk -F'\t' -v OFS='\t' -v n=${SMP1_START} '{ $2 -= n; $3 -= n }1' ${OUTDIR}/${SMP_PAIR}_BED2.bed9 | grep ${SMP1} > ${LOCAL_FOLDER}/${SMP1}_sv_horhap.bed

awk -F'\t' -v OFS='\t' -v n=${SMP2_START} '{ $2 -= n; $3 -= n }1' ${OUTDIR}/${SMP_PAIR}_BED2.bed9 | grep ${SMP2} > ${LOCAL_FOLDER}/${SMP2}_sv_horhap.bed

# CENTROLIGN
awk -F'\t' -v OFS='\t' -v n=${SMP1_START} '{ $2 -= n; $3 -= n }1' ${OUTDIR}/${SMP_PAIR}_BED1.bed9 | grep ${SMP1} > ${LOCAL_FOLDER}/${SMP1}_sv_cen.bed

awk -F'\t' -v OFS='\t' -v n=${SMP2_START} '{ $2 -= n; $3 -= n }1' ${OUTDIR}/${SMP_PAIR}_BED1.bed9 | grep ${SMP2} > ${LOCAL_FOLDER}/${SMP2}_sv_cen.bed

# Generate synteny plots
source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate synteny

mkdir -p ${OUTDIR}/synteny_plots
# horhap view
python /private/groups/migalab/juklucas/centrolign/chr12_test125/synteny_plot_bokeh.py   \
    --beds \
        ${LOCAL_FOLDER}/${SMP1}_sv_horhap.bed \
        ${LOCAL_FOLDER}/${SMP2}_sv_horhap.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrX/induced_pairwise_cigars/pairwise_cigar_${SMP1}_${SMP2}.txt \
    --output ${OUTDIR}/synteny_plots/md3k_msf0.5_${SMP1}_${SMP2}_synteny.horhap.html \
    --show-mismatches \
    --web

# centrolign view
python /private/groups/migalab/juklucas/centrolign/chr12_test125/synteny_plot_bokeh.py   \
    --beds \
        ${LOCAL_FOLDER}/${SMP1}_sv_cen.bed \
        ${LOCAL_FOLDER}/${SMP2}_sv_cen.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrX/induced_pairwise_cigars/pairwise_cigar_${SMP1}_${SMP2}.txt \
    --output ${OUTDIR}/synteny_plots/md3k_msf0.5_${SMP1}_${SMP2}_synteny.cen.html \
    --show-mismatches \
    --web

rm -rf ${LOCAL_FOLDER}
```

Submitting script
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/f1_scores/chrX/sv_compare2
cd /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/f1_scores/chrX/sv_compare2

mkdir md3k_msf0.5
mkdir -p logs

sbatch sv_compare.sh
```

## Trying chrX with 2x average STV size (4500)
```sh
#!/bin/bash
#SBATCH --job-name=sv_compare2_chrX
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --output=logs/sv_compare_%x.%j.log
#SBATCH --array=[1-6991]%100

BED_CSV=/private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/f1_scores/chrX/matched_beds.filt.csv

CENTROLIGN_BED=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$BED_CSV" | cut -f1 -d",")
HORHAP_BED=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$BED_CSV" | cut -f2 -d",")

SMP_PAIR=$(basename -s .bed $CENTROLIGN_BED )
SMP1="${SMP_PAIR%%_*}"
SMP2="${SMP_PAIR##*_}"  

echo $SMP_PAIR
echo $CENTROLIGN_BED
echo $HORHAP_BED

LOCAL_FOLDER=/data/tmp/$(whoami)/chrX_${SMP_PAIR}/
mkdir -p ${LOCAL_FOLDER}

OUTDIR=/private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/f1_scores/chrX/sv_compare2/md4.5k_msf0.5
mkdir -p $OUTDIR/

# select only insertions from query, deletions from ref coords
grep "I" ${CENTROLIGN_BED} | cut -f4,5,6,8 > ${LOCAL_FOLDER}/${SMP_PAIR}.cen.bed

grep "D" ${CENTROLIGN_BED} | cut -f1-3,8 >> ${LOCAL_FOLDER}/${SMP_PAIR}.cen.bed

grep "I" ${HORHAP_BED} | cut -f4,5,6,8 > ${LOCAL_FOLDER}/${SMP_PAIR}.hor.bed

grep "D" ${HORHAP_BED} | cut -f1-3,8 >> ${LOCAL_FOLDER}/${SMP_PAIR}.hor.bed

# Run SV comparison
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/sv_comparev2.py \
    ${LOCAL_FOLDER}/${SMP_PAIR}.cen.bed \
    ${LOCAL_FOLDER}/${SMP_PAIR}.hor.bed \
    --max_dist 4500 --min_size_frac 0.5 \
    --out_prefix ${OUTDIR}/${SMP_PAIR} --bed9 --debug

# Auto-generate synteny plots showing cenhap and horhap SV qualities
# get start coord of alpha sat sequence
SMP1_START=`grep chrX /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/${SMP1}_asat_arrays.bed | cut -f2`

SMP2_START=`grep chrX /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/${SMP2}_asat_arrays.bed | cut -f2`

echo $SMP1_START
echo $SMP2_START

# HORHAP
awk -F'\t' -v OFS='\t' -v n=${SMP1_START} '{ $2 -= n; $3 -= n }1' ${OUTDIR}/${SMP_PAIR}_BED2.bed9 | grep ${SMP1} > ${LOCAL_FOLDER}/${SMP1}_sv_horhap.bed

awk -F'\t' -v OFS='\t' -v n=${SMP2_START} '{ $2 -= n; $3 -= n }1' ${OUTDIR}/${SMP_PAIR}_BED2.bed9 | grep ${SMP2} > ${LOCAL_FOLDER}/${SMP2}_sv_horhap.bed

# CENTROLIGN
awk -F'\t' -v OFS='\t' -v n=${SMP1_START} '{ $2 -= n; $3 -= n }1' ${OUTDIR}/${SMP_PAIR}_BED1.bed9 | grep ${SMP1} > ${LOCAL_FOLDER}/${SMP1}_sv_cen.bed

awk -F'\t' -v OFS='\t' -v n=${SMP2_START} '{ $2 -= n; $3 -= n }1' ${OUTDIR}/${SMP_PAIR}_BED1.bed9 | grep ${SMP2} > ${LOCAL_FOLDER}/${SMP2}_sv_cen.bed

# Generate synteny plots
source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate synteny

mkdir -p ${OUTDIR}/synteny_plots
# horhap view
python /private/groups/migalab/juklucas/centrolign/chr12_test125/synteny_plot_bokeh.py   \
    --beds \
        ${LOCAL_FOLDER}/${SMP1}_sv_horhap.bed \
        ${LOCAL_FOLDER}/${SMP2}_sv_horhap.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrX/induced_pairwise_cigars/pairwise_cigar_${SMP1}_${SMP2}.txt \
    --output ${OUTDIR}/synteny_plots/md4.5k_msf0.5_${SMP1}_${SMP2}_synteny.horhap.html \
    --show-mismatches \
    --web

# centrolign view
python /private/groups/migalab/juklucas/centrolign/chr12_test125/synteny_plot_bokeh.py   \
    --beds \
        ${LOCAL_FOLDER}/${SMP1}_sv_cen.bed \
        ${LOCAL_FOLDER}/${SMP2}_sv_cen.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrX/induced_pairwise_cigars/pairwise_cigar_${SMP1}_${SMP2}.txt \
    --output ${OUTDIR}/synteny_plots/md4.5k_msf0.5_${SMP1}_${SMP2}_synteny.cen.html \
    --show-mismatches \
    --web

rm -rf ${LOCAL_FOLDER}
```

Submitting script
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/f1_scores/chrX/sv_compare2
cd /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/f1_scores/chrX/sv_compare2

mkdir md4.5k_msf0.5
mkdir -p logs

sbatch sv_compare.sh
```
