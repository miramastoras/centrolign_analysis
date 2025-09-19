### Spot checking some alignments for chr3 and chr4, which have an embedded HSAT array

Select some samples to check
```
# distance matrix files
cd /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_all_pairs
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/distance_matrices

HG00642.2,HG01175.1,0.5073888442241132
HG01175.1,NA19036.1,0.06851893599718117
NA19036.1,NA19159.2,0.5247622937758846
NA19159.2,NA21110.2,0.736252699915496
NA21110.2,NA21144.1,0.9561056655670755
NA21144.1,NA21144.2,0.6491349614679709

```
Run some dotplots:
```

time docker run \
  -u `id -u`:`id -g` \
  -v "/private/groups/patenlab/mira":"/private/groups/patenlab/mira" \
  miramastoras/centromere_scripts:v0.1.2 \
  python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/visualization/plot_dotplot_alignment.py \
  /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/HG00128_hap2/analysis/extract_hors_HPRC_outputs/HG00128_hap2_HG00128.2_chr3_hor_array.fasta \
  /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/NA19036_hap1/analysis/extract_hors_HPRC_outputs/NA19036_hap1_NA19036.1_chr3_hor_array.fasta \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr3/pairwise_cigar/pairwise_cigar_HG00128.2_NA19036.1.txt \
  /private/groups/patenlab/mira/centrolign/analysis/chr3_chr4_interspersed_hsat_r2/HG00128.2_NA19036.1_chr3.release2.dotplot.svg

s
time docker run \
  -u `id -u`:`id -g` \
  -v "/private/groups/patenlab/mira":"/private/groups/patenlab/mira" \
  miramastoras/centromere_scripts:v0.1.2 \
  python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/visualization/plot_dotplot_alignment.py \
  /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/NA19036_hap1/analysis/extract_hors_HPRC_outputs/NA19036_hap1_NA19036.1_chr3_hor_array.fasta \
  /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/NA19159_hap2/analysis/extract_hors_HPRC_outputs/NA19159_hap2_NA19159.2_chr3_hor_array.fasta \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr3/pairwise_cigar/pairwise_cigar_NA19036.1_NA19159.2.txt \
  /private/groups/patenlab/mira/centrolign/analysis/chr3_chr4_interspersed_hsat_r2/NA19036.1_NA19159.2_chr3.release2.dotplot.svg

/private/groups/patenlab/mira/centrolign/analysis/chr3_chr4_interspersed_hsat_r2

time docker run \
  -u `id -u`:`id -g` \
  -v "/private/groups/patenlab/mira":"/private/groups/patenlab/mira" \
  miramastoras/centromere_scripts:v0.1.2 \
  python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/visualization/plot_dotplot_alignment.py \
  /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/NA19036_hap1/analysis/extract_hors_HPRC_outputs/NA19036_hap1_NA19036.1_chr3_hor_array.fasta \
  /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/NA19159_hap2/analysis/extract_hors_HPRC_outputs/NA19159_hap2_NA19159.2_chr3_hor_array.fasta \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr3/pairwise_cigar/pairwise_cigar_NA19036.1_NA19159.2.txt \
  /private/groups/patenlab/mira/centrolign/analysis/chr3_chr4_interspersed_hsat_r2/NA19036.1_NA19159.2_chr3.release2.dotplot.svg


HG01993.2,HG02273.1,0.19277374877997555
time docker run \
  -u `id -u`:`id -g` \
  -v "/private/groups/patenlab/mira":"/private/groups/patenlab/mira" \
  miramastoras/centromere_scripts:v0.1.2 \
  python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/visualization/plot_dotplot_alignment.py \
  /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/HG01993_mat/analysis/extract_hors_HPRC_outputs/HG01993_mat_HG01993.2_chr3_hor_array.fasta \
  /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/HG02273_hap1/analysis/extract_hors_HPRC_outputs/HG02273_hap1_HG02273.1_chr3_hor_array.fasta \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr3/pairwise_cigar/pairwise_cigar_HG01993.2_HG02273.1.txt \
  /private/groups/patenlab/mira/centrolign/analysis/chr3_chr4_interspersed_hsat_r2/HG01993.2_HG02273.1_chr3.release2.dotplot.svg


HG00126.2,HG03239.1,0.4850908675003496
time docker run \
  -u `id -u`:`id -g` \
  -v "/private/groups/patenlab/mira":"/private/groups/patenlab/mira" \
  miramastoras/centromere_scripts:v0.1.2 \
  python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/visualization/plot_dotplot_alignment.py \
  /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/HG00126_hap2/analysis/extract_hors_HPRC_outputs/HG00126_hap2_HG00126.2_chr4_hor_array.fasta \
  /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/HG03239_pat/analysis/extract_hors_HPRC_outputs/HG03239_pat_HG03239.1_chr4_hor_array.fasta \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr4/pairwise_cigar/pairwise_cigar_HG00126.2_HG03239.1.txt \
  /private/groups/patenlab/mira/centrolign/analysis/chr3_chr4_interspersed_hsat_r2/HG00126.2_HG03239.1_chr3.release2.dotplot.svg

```

Run Julian's alignment visualization - chr3
```
cd /private/groups/patenlab/mira/centrolign/analysis/chr3_chr4_interspersed_hsat_r2

conda activate aws

sed -i 's/\r$//' bed_files_synteny.csv

# intersect censat bed file with centrolign bed file, reset to 0
while IFS=',' read -r smp col2 col3 censatBed centrolignBed; do
    echo $smp
    aws s3 cp $censatBed censat_beds/
    censat_base=`basename ${censatBed}`

    grep chr3 ${centrolignBed} | bedtools intersect -a censat_beds/${censat_base} -b - | awk 'NR==1 {offset=$2} {print $1, $2 - offset, $3 - offset, $4, $5, $6, $7, $8, $9}' OFS='\t' > censat_beds/${smp}_cenSat.reset.bed

    done < bed_files_synteny.csv
    | awk 'NR==1 {offset=$2} {print $1, $2 - offset, $3 - offset, $4, $5, $6, $7, $8, $9}' OFS='\t' > censat_beds/${smp}_cenSat.reset.bed
done < bed_files_synteny.csv

python /private/groups/migalab/juklucas/centrolign/chr12_test125/synteny_plot_bokeh.py   \
    --beds \
      censat_beds/HG00642_mat_cenSat.reset.bed \
      censat_beds/HG01175_pat_cenSat.reset.bed \
      censat_beds/NA19036_hap1_cenSat.reset.bed \
      censat_beds/NA19159_hap2_cenSat.reset.bed \
      censat_beds/NA21110_hap2_cenSat.reset.bed \
      censat_beds/NA21144_hap1_cenSat.reset.bed \
      censat_beds/NA21144_hap2_cenSat.reset.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr3/pairwise_cigar/pairwise_cigar_HG00642.2_HG01175.1.txt \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr3/pairwise_cigar/pairwise_cigar_HG01175.1_NA19036.1.txt \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr3/pairwise_cigar/pairwise_cigar_NA19036.1_NA19159.2.txt \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr3/pairwise_cigar/pairwise_cigar_NA19159.2_NA21110.2.txt \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr3/pairwise_cigar/pairwise_cigar_NA21110.2_NA21144.1.txt \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr3/pairwise_cigar/pairwise_cigar_NA21144.1_NA21144.2.txt \
    --output synteny_by_identity_chr3_rel2.html \
    --web
```

Run synteny plots for Chr 4:

Selecting samples:
```
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/distance_matrices

HG00099.2,HG00128.1,0.5692303284950344
HG00128.1,HG00146.2,0.38878271263813824
HG00146.2,HG01530.2,0.710622275131364
HG01530.2,NA21102.2,0.693915995504394
```

```
cd /private/groups/patenlab/mira/centrolign/analysis/chr3_chr4_interspersed_hsat_r2


while IFS=',' read -r smp col2 col3 censatBed centrolignBed; do
    echo $smp
    aws s3 cp $censatBed censat_beds/
    censat_base=`basename ${censatBed}`

    grep chr4 ${centrolignBed} | bedtools intersect -a censat_beds/${censat_base} -b - | awk 'NR==1 {offset=$2} {print $1, $2 - offset, $3 - offset, $4, $5, $6, $7, $8, $9}' OFS='\t' > censat_beds/${smp}_cenSat.reset.chr4.bed

  done < bedfiles_synteny_chr4.csv


python /private/groups/migalab/juklucas/centrolign/chr12_test125/synteny_plot_bokeh.py   \
    --beds \
      censat_beds/HG00099_hap2_cenSat.reset.chr4.bed \
      censat_beds/HG00128_hap1_cenSat.reset.chr4.bed \
      censat_beds/HG00146_hap2_cenSat.reset.chr4.bed \
      censat_beds/HG01530_hap2_cenSat.reset.chr4.bed \
      censat_beds/NA21102_hap2_cenSat.reset.chr4.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr4/pairwise_cigar/pairwise_cigar_HG00099.2_HG00128.1.txt \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr4/pairwise_cigar/pairwise_cigar_HG00128.1_HG00146.2.txt \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr4/pairwise_cigar/pairwise_cigar_HG00146.2_HG01530.2.txt \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr4/pairwise_cigar/pairwise_cigar_HG01530.2_NA21102.2.txt \
    --output synteny_by_identity_chr4_rel2.html \
    --web
```
