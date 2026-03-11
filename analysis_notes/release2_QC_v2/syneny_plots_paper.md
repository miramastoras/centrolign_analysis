### Code for final synteny plots for centrolign censat_paper

sample pairs from HSAT 3 and 4 with dist < 0.2
```sh
/private/groups/patenlab/mira/centrolign/analysis/chr3_chr4_interspersed_hsat_r2/

chr3
HG01928.1,HG01978.2
HG01978.2,HG01993.2
HG01993.2,HG02293.2
HG02293.2,NA19776.2

cd /private/groups/patenlab/mira/centrolign/analysis/chr3_chr4_interspersed_hsat_r2
source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate base


SMPS=("HG01928.1" "HG01978.2" "HG01993.2" "HG02293.2" "NA19776.2")

for smp in "${SMPS[@]}"
do
  ID=`echo $smp | sed 's/\./,/g'`

  BED=`grep $ID /private/groups/patenlab/mira/centrolign/annotations/censat_hprc_r2_v1.0.index.csv | cut -f4 -d","`
  echo $BED

  aws s3 cp $BED censat_beds_r2_QC_v2/
  echo $censat_base
  censat_base=`basename ${BED}`

  grep chr3 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/${smp}_asat_arrays.bed | bedtools intersect -a censat_beds_r2_QC_v2/${censat_base} -b - | awk 'NR==1 {offset=$2} {print $1, $2 - offset, $3 - offset, $4, $5, $6, $7, $8, $9}' OFS='\t' > censat_beds_r2_QC_v2/${smp}_cenSat.reset.chr3.bed

done
```

Figure panel ready version with legends:
```sh
python /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/synteny_plot_static.py  \
    --beds \
        /private/groups/patenlab/mira/centrolign/analysis/chr3_chr4_interspersed_hsat_r2/censat_beds_r2_QC_v2/HG01928.1_cenSat.reset.chr3.bed \
        /private/groups/patenlab/mira/centrolign/analysis/chr3_chr4_interspersed_hsat_r2/censat_beds_r2_QC_v2/HG01978.2_cenSat.reset.chr3.bed \
        /private/groups/patenlab/mira/centrolign/analysis/chr3_chr4_interspersed_hsat_r2/censat_beds_r2_QC_v2/HG01993.2_cenSat.reset.chr3.bed \
        /private/groups/patenlab/mira/centrolign/analysis/chr3_chr4_interspersed_hsat_r2/censat_beds_r2_QC_v2/HG02293.2_cenSat.reset.chr3.bed \
        /private/groups/patenlab/mira/centrolign/analysis/chr3_chr4_interspersed_hsat_r2/censat_beds_r2_QC_v2/NA19776.2_cenSat.reset.chr3.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr3/induced_pairwise_cigars/pairwise_cigar_HG01928.1_HG01978.2.txt \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr3/induced_pairwise_cigars/pairwise_cigar_HG01978.2_HG01993.2.txt \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr3/induced_pairwise_cigars/pairwise_cigar_HG01993.2_HG02293.2.txt \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr3/induced_pairwise_cigars/pairwise_cigar_HG02293.2_NA19776.2.txt \
    --show-mismatches \
    --labels "HG01928 hap1" "HG01978 hap2" "HG01993 hap2" "HG02293 hap2" "NA19776 hap2" \
    --output /private/groups/patenlab/mira/centrolign/analysis/chr3_chr4_interspersed_hsat_r2/release2_QC_v2_figs/Chr3_embeddedHsat_synteny.svg
```

### HOR HAP sv concordance examples

### Main panel synteny plot

```sh
conda activate synteny

cut -f4-6 /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr11/SV_beds/chr11_subgroupA/HG01243.1_HG01884.1.bed > /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01884.1_HG01243.1.centrolign_SVs.bed

# File 1 (8 columns) — use col 1-3 as chrom/start/end, add dummy BED9 cols
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, ".", 0, ".", $2, $3, "0,0,0"}' \
    /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/chr11/SV_beds/chr11_subgroupA/HG01243.1_HG01884.1.bed \
    > /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01243.1_HG01884.1.centrolign_SVs.bed9

# File 2 (3 columns) — same
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, ".", 0, ".", $2, $3, "0,0,0"}' \
    /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01884.1_HG01243.1.centrolign_SVs.bed \
    > /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01884.1_HG01243.1.centrolign_SVs.bed9


python /private/groups/migalab/juklucas/centrolign/chr12_test125/synteny_plot_bokeh.py   \
    --beds \
        /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01243.1_HG01884.1.centrolign_SVs.bed9 \
        /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01884.1_HG01243.1.centrolign_SVs.bed9 \
    --cigars \
         \
    --output /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01243.1_HG01884.1.centrolign.synteny.html \
    --show-mismatches \
    --web

# convert horhap sv bed to array coords
grep chr11 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/HG01243.1_asat_arrays.bed
# HG01243#1#CM099608.1	50804872	53964932	chr11

awk 'BEGIN{OFS="\t"} {$2=($2-50804872<0?0:$2-50804872); $3=($3-50804872<0?0:$3-50804872); print}' /private/groups/migalab/fryabov/HPRC/horhap_alignments/chr11/CV_bed/HG01243.1_HG01884.1.bed | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, ".", 0, ".", $2, $3, "0,0,0"}' > /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01243.1_HG01884.1.array_coords.horhap.bed

grep chr11 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/HG01884.1_asat_arrays.bed
# HG01884#1#CM086410.1	50951635	54608852	chr11 # 3657217

cut -f 4-6 /private/groups/migalab/fryabov/HPRC/horhap_alignments/chr11/CV_bed/HG01243.1_HG01884.1.bed | awk 'BEGIN{OFS="\t"} {$2=($2-50951635<0?0:$2-50951635); $3=($3-50951635<0?0:$3-50951635); print}' | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, ".", 0, ".", $2, $3, "0,0,0"}' > /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01884.1_HG01243.1.array_coords.horhap.bed


# convert M into = for cigar
sed 's/M/=/g' /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01243.1_HG01884.1.horhap_cigar.txt > /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01243.1_HG01884.1.horhap_cigar_converted.txt
```

Horhap beds:
```sh
awk 'BEGIN{OFS="\t"} {$2=($2-50804872<0?0:$2-50804872); $3=($3-50804872<0?0:$3-50804872); print}' /private/groups/migalab/fryabov/HPRC/horhaps/data/chr11/k5/HG01243#1#CM099608.1.bed > /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01243.1_horhaps.array_coords.bed

awk 'BEGIN{OFS="\t"} {$2=($2-50951635<0?0:$2-50951635); $3=($3-50951635<0?0:$3-50951635); print}' /private/groups/migalab/fryabov/HPRC/horhaps/data/chr11/k5/HG01884#1#CM086410.1.bed > /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01884.1_horhaps.array_coords.bed
```

Horhap colors, flipped orientation

```sh
python3 /private/groups/patenlab/mira/centrolign/github/censat_paper/scripts/centrolign_result_parsing/reverse_cigar.py         /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01243.1_HG01884.1.horhap_cigar.txt


python /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/synteny_plot_static.py  \
    --beds \
        /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01884.1_horhaps.array_coords.bed \
        /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01243.1_horhaps.array_coords.bed \
        /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01884.1_horhaps.array_coords.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01243.1_HG01884.1.horhap_cigar_reversed.txt \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr11/subgroup_A/induced_pairwise_cigars/pairwise_cigar_HG01243.1_HG01884.1.txt \
    --labels "HG01884 hap1" "HG01243 hap1" "HG01884 hap1" \
    --alignment-labels "HORHap align" "Centrolign" \
    --output /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01243.1_HG01884.1.final.main_panel.png
```

### Supplementary figures synteny plot


#### Worse example for chr11

HG01433.2 HG01346.1

Get asat start coords
```sh
# convert horhap sv bed to array coords
grep chr11 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/HG01433.2_asat_arrays.bed
# HG01433#2#CM086516.1	50858424	54265807	chr11

grep chr11 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/HG01346.1_asat_arrays.bed
# HG01346#1#CM086591.1	50697180	54102506	chr11
```

Convert horhap beds to array coords
```sh
awk 'BEGIN{OFS="\t"} {$2=($2-50858424<0?0:$2-50858424); $3=($3-50858424<0?0:$3-50858424); print}' /private/groups/migalab/fryabov/HPRC/horhaps/data/chr11/k5/HG01433#2#CM086516.1.bed > /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01433.2_horhaps.array_coords.bed

awk 'BEGIN{OFS="\t"} {$2=($2-50697180<0?0:$2-50697180); $3=($3-50697180<0?0:$3-50697180); print}' /private/groups/migalab/fryabov/HPRC/horhaps/data/chr11/k5/HG01346#1#CM086591.1.bed > /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01346.1_horhaps.array_coords.bed
```

Get horhap cigar in right format
```sh
grep HG01433.2 /private/groups/migalab/fryabov/HPRC/horhap_alignments/chr11/cigar_strings.tsv | grep HG01346.1 | cut -f 2 > /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01433.2_HG01346.1.horhap_cigar.txt

python3 /private/groups/patenlab/mira/centrolign/github/censat_paper/scripts/centrolign_result_parsing/reverse_cigar.py /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01433.2_HG01346.1.horhap_cigar.txt

sed 's/M/=/g' /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01433.2_HG01346.1.horhap_cigar_reversed.txt > /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01433.2_HG01346.1.horhap_cigar_reversed_converted.txt
```

Synteny
```sh
python /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/synteny_plot_static.py \
    --beds \
        /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01346.1_horhaps.array_coords.bed \
        /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01433.2_horhaps.array_coords.bed \
        /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01346.1_horhaps.array_coords.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01433.2_HG01346.1.horhap_cigar_reversed_converted.txt \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr11/subgroup_B/induced_pairwise_cigars/pairwise_cigar_HG01433.2_HG01346.1.txt \
    --labels "HG01346 hap1" "HG01433 hap2" "HG01346 hap1" \
    --colorpop \
    --alignment-labels "HORHap align" "Centrolign" \
    --output /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01433.2_HG01346.1.supp.final.svg
```

##### Example 2

HG01975.1 HG01928.2 0.189737   0.701031 0.576037 0.632417


Get asat start coords
```sh
# convert horhap sv bed to array coords
grep chr11 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/HG01975.1_asat_arrays.bed
# HG01975#1#CM099802.1	50938159	53773835	chr11

grep chr11 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/HG01928.2_asat_arrays.bed
# HG01928#2#CM089009.1	50953069	53619097	chr11
```

Convert horhap beds to array coords
```sh
awk 'BEGIN{OFS="\t"} {$2=($2-50938159<0?0:$2-50938159); $3=($3-50938159<0?0:$3-50938159); print}' /private/groups/migalab/fryabov/HPRC/horhaps/data/chr11/k5/HG01975#1#CM099802.1.bed > /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01975.1_horhaps.array_coords.bed

awk 'BEGIN{OFS="\t"} {$2=($2-50953069<0?0:$2-50953069); $3=($3-50953069<0?0:$3-50953069); print}' /private/groups/migalab/fryabov/HPRC/horhaps/data/chr11/k5/HG01928#2#CM089009.1.bed > /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01928.2_horhaps.array_coords.bed
```

Get horhap cigar in right format
```sh
grep HG01928.2 /private/groups/migalab/fryabov/HPRC/horhap_alignments/chr11/cigar_strings.tsv | grep HG01975.1 | cut -f 2 > /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01975.1_HG01928.2.horhap_cigar.txt

python3 /private/groups/patenlab/mira/centrolign/github/censat_paper/scripts/centrolign_result_parsing/reverse_cigar.py  /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01975.1_HG01928.2.horhap_cigar.txt

sed 's/M/=/g' /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01975.1_HG01928.2.horhap_cigar_reversed.txt > /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01975.1_HG01928.2.horhap_cigar_reversed_converted.txt
```

Synteny
```sh
python /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/synteny_plot_static.py  \
    --beds \
        /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01928.2_horhaps.array_coords.bed \
        /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01975.1_horhaps.array_coords.bed \
        /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01928.2_horhaps.array_coords.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01975.1_HG01928.2.horhap_cigar_reversed_converted.txt \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr11/subgroup_A/induced_pairwise_cigars/pairwise_cigar_HG01975.1_HG01928.2.txt \
    --labels "HG01928 hap2" "HG01975 hap1" "HG01928 hap2" \
    --colorpop \
    --alignment-labels "HORHap align" "Centrolign" \
    --output /private/groups/patenlab/mira/centrolign/analysis/horhap_SV_concordance/synteny_figure/HG01975.1_HG01928.2.sup.final.svg
```
