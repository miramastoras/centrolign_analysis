## Investigating nucleotide diversity and Ti/Tv ratio for HPRC r2 censat alpha arrays

### Running Jordan's variant matrix script and Nucleotide diversity / Ti/Tv ratio with default settings

Prepare input file for running slurm array
```sh
cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/centrolign_results.csv | grep "release 2 QC v2" | grep -v "chr20," | cut -f1,9,10 -d","   | while IFS=',' read -r subgroup gfa induced ; do
  echo $subgroup,$gfa,$induced
done > /private/groups/patenlab/mira/centrolign/analysis/variant_matrix/make_var_mat_input_index.csv
```

```sh
#!/bin/bash
#SBATCH --job-name=make_var_mat
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x.%j.log
#SBATCH --time=12:00:00
#SBATCH --array=[1]%43

IDX_FILE="/private/groups/patenlab/mira/centrolign/analysis/variant_matrix/make_var_mat_input_index.csv"
OUT_BASE="/private/groups/patenlab/mira/centrolign/analysis/variant_matrix/default"

CHR_SUBGROUP=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$IDX_FILE" | cut -f1 -d",")
GFA=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$IDX_FILE" | cut -f2 -d",")
INDUCED_PATH=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$IDX_FILE" | cut -f3 -d",")
OUTDIR=${OUT_BASE}/${CHR_SUBGROUP}

echo $IDX_FILE
echo $OUT_BASE
echo $CHR_SUBGROUP
echo $GFA
echo $INDUCED_PATH
echo $OUTDIR

mkdir -p $OUTDIR

~/progs/centrolign/build/make_var_mat -b \
    $GFA \
    > ${OUTDIR}/${CHR_SUBGROUP}.var_mat.tsv

python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/estimate_nuc_diversity.py \
    ${OUTDIR}/${CHR_SUBGROUP}.var_mat.tsv \
    ${INDUCED_PATH}/pairwise_cigar_ \
    &> ${OUTDIR}/${CHR_SUBGROUP}.nucdiv_results.txt
```

### Subsetting calculation to pairs < 0.2 and < 0.1

Implementation and unit testing
```sh
# implemented in here
/private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/estimate_nuc_diversity_v2.py

# Claude unit tests
/private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/test_estimate_nuc_diversity_v2.py
/private/groups/patenlab/mira/centrolign/analysis/variant_matrix/unit_tests
```

```sh
#!/bin/bash
#SBATCH --job-name=test_nuc_div_update
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x.%j.log
#SBATCH --time=12:00:00

# Testing with a pairwise dist threshold of 1 to confirm matching results with prev implementation
time python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/estimate_nuc_diversity_v2.py \
    /private/groups/patenlab/mira/centrolign/analysis/variant_matrix/default/chr12_subgroup1/chr12_subgroup1.var_mat.tsv \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/induced_pairwise_cigars/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/chr12_r2_QC_v2_centrolign_pairwise_distance.csv \
    --threshold 1.2 \
    &> /private/groups/patenlab/mira/centrolign/analysis/variant_matrix/unit_tests/chr12_subgroup1.nucdiv_results.txt

#!/bin/sh
# Testing with a pairwise dist threshold of 1 to confirm matching results with prev implementation
time python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/estimate_nuc_diversity_v2.py \
    /private/groups/patenlab/mira/centrolign/analysis/variant_matrix/default/chr14_subgroup_A/chr14_subgroup_A.var_mat.tsv \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr14/subgroup_A/induced_pairwise_cigars/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/chr14_r2_QC_v2_centrolign_pairwise_distance.csv \
    --threshold 1.2 \
    &> /private/groups/patenlab/mira/centrolign/analysis/variant_matrix/unit_tests/chr14_subgroupA.nucdiv_results.txt
```

Running on pairs less than 0.1 and 0.2
```sh
#!/bin/bash
#SBATCH --job-name=make_var_mat_dist0.2
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x.%j.log
#SBATCH --time=1:00:00
#SBATCH --array=[1-60]%60

IDX_FILE="/private/groups/patenlab/mira/centrolign/analysis/variant_matrix/make_var_mat_input_index.csv"
OUT_BASE="/private/groups/patenlab/mira/centrolign/analysis/variant_matrix/pairs_lt_0.2"
MAT_PATH="/private/groups/patenlab/mira/centrolign/analysis/variant_matrix/default"

CHR_SUBGROUP=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$IDX_FILE" | cut -f1 -d",")
CHR=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$IDX_FILE" | cut -f1 -d"," | cut -f1 -d"_")
GFA=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$IDX_FILE" | cut -f2 -d",")
INDUCED_PATH=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$IDX_FILE" | cut -f3 -d",")
OUTDIR=${OUT_BASE}/${CHR_SUBGROUP}

echo $IDX_FILE
echo $OUT_BASE
echo $CHR_SUBGROUP
echo $CHR
echo $GFA
echo $INDUCED_PATH
echo $OUTDIR

mkdir -p $OUTDIR

python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/estimate_nuc_diversity_v2.py \
    ${MAT_PATH}/${CHR_SUBGROUP}/${CHR_SUBGROUP}.var_mat.tsv \
    ${INDUCED_PATH}/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/${CHR}_r2_QC_v2_centrolign_pairwise_distance.csv \
    --threshold 0.2 \
    &> ${OUTDIR}/${CHR_SUBGROUP}.nucdiv_results.txt
```

```sh
#!/bin/bash
#SBATCH --job-name=make_var_mat_dist0.1
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x.%j.log
#SBATCH --time=1:00:00
#SBATCH --array=[1-60]%60

IDX_FILE="/private/groups/patenlab/mira/centrolign/analysis/variant_matrix/make_var_mat_input_index.csv"
OUT_BASE="/private/groups/patenlab/mira/centrolign/analysis/variant_matrix/pairs_lt_0.1"
MAT_PATH="/private/groups/patenlab/mira/centrolign/analysis/variant_matrix/default"

CHR_SUBGROUP=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$IDX_FILE" | cut -f1 -d",")
CHR=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$IDX_FILE" | cut -f1 -d"," | cut -f1 -d"_")
GFA=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$IDX_FILE" | cut -f2 -d",")
INDUCED_PATH=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$IDX_FILE" | cut -f3 -d",")
OUTDIR=${OUT_BASE}/${CHR_SUBGROUP}

echo $IDX_FILE
echo $OUT_BASE
echo $CHR_SUBGROUP
echo $CHR
echo $GFA
echo $INDUCED_PATH
echo $OUTDIR

mkdir -p $OUTDIR

python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/estimate_nuc_diversity_v2.py \
    ${MAT_PATH}/${CHR_SUBGROUP}/${CHR_SUBGROUP}.var_mat.tsv \
    ${INDUCED_PATH}/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/${CHR}_r2_QC_v2_centrolign_pairwise_distance.csv \
    --threshold 0.1 \
    &> ${OUTDIR}/${CHR_SUBGROUP}.nucdiv_results.txt
```
### New version of script to filter any SNVs within 20 bp of an indel

```sh
~/progs/centrolign/build/make_var_mat -b \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr5/subgroup_0/chr5.subgroup_0.centrolign.gfa \
    --svs --indels --allow-nest --mnvs --positions \
    > /private/groups/patenlab/mira/centrolign/analysis/variant_matrix/unit_tests/chr5_subgroup_0.svs_indels_mnvs.var_mat.tsv


python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/estimate_nuc_diversity_v3.py \
    /private/groups/patenlab/mira/centrolign/analysis/variant_matrix/unit_tests/chr5_subgroup_0.svs_indels_mnvs.var_mat.tsv \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr5/subgroup_0/induced_pairwise_cigars/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/chr5_r2_QC_v2_centrolign_pairwise_distance.csv \
    --threshold 1.2 \
    --indel_window 0
```


### Run pairwise alignments using Unialigner, Rama, and Minimap2 on sample pairs < 0.2, to orthogonally validate the low Ti/Tv ratio we are seeing

Create index file of all the distance matrices with pairwise alignments < 0.2
```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices_lt0.2

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")
for chr in "${chromosomes[@]}"
do
    awk -v chr="$chr" 'BEGIN{OFS=","} {print $0, chr}' ${chr}_r2_QC_v2_centrolign_pairwise_distance.csv
done >> /private/groups/patenlab/mira/centrolign/analysis/variant_matrix/all_pairs_lt0.2_index.csv
```

```sh
#!/bin/bash
#SBATCH --job-name=orthogonal_aligners
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --output=logs/%x.%j.log
#SBATCH --time=1:00:00
#SBATCH --array=[1]%150

# Parse row in index file
IDX_FILE="/private/groups/patenlab/mira/centrolign/analysis/variant_matrix/all_pairs_lt0.2_index.csv"
SAMP1=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$IDX_FILE" | cut -f1 -d",")
SAMP2=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$IDX_FILE" | cut -f2 -d",")
DIST=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$IDX_FILE" | cut -f3 -d",")
CHR=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$IDX_FILE" | cut -f4 -d",")

# Define software paths
UNIALIGNER=/private/groups/patenlab/mira/centrolign/github/unialigner/tandem_aligner/build/bin/tandem_aligner

# Define fasta paths
FASTA1=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/${CHR}/${SAMP1}_${CHR}_hor_array.fasta
FASTA2=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/${CHR}/${SAMP2}_${CHR}_hor_array.fasta

echo "aligning with unaligner"
UNI_OUTFILE=/private/groups/patenlab/mira/centrolign/analysis/variant_matrix/orthogonal_aligners/unialigner/${CHR}/unialigner_cigar_${SAMP1}_${SAMP2}.txt
mkdir -p `dirname $UNI_OUTFILE`
UNI_OUT_TMP=/private/groups/patenlab/mira/centrolign/analysis/variant_matrix/orthogonal_aligners/unialigner/${CHR}/${SAMP1}_${SAMP2}_tmp/

ulimit -v 471859200 && /usr/bin/time -v ${UNIALIGNER} --first $FASTA1 --second $FASTA2 -o $UNI_OUT_TMP
mv $UNI_OUT_TMP/cigar.txt $UNI_OUTFILE
# delete the rest of the output
rm -r $UNI_OUT_TMP

echo "aligning with RAMA"
RAMA_TEMP_OUTDIR=/private/groups/patenlab/mira/centrolign/analysis/variant_matrix/orthogonal_aligners/Rama/${CHR}/${SAMP1}_${SAMP2}_tmp/
RAMA_OUTFILE=/private/groups/patenlab/mira/centrolign/analysis/variant_matrix/orthogonal_aligners/Rama/${CHR}/Rama_cigar_${SAMP1}_${SAMP2}.txt
mkdir -p `dirname $RAMA_OUTFILE`

ulimit -v 471859200 && /usr/bin/time -v \
docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
    miramastoras/rama:latest ./RaMA -t 5 \
    -r $FASTA1 \
    -q $FASTA2 \
    -o $RAMA_TEMP_OUTDIR

mv $RAMA_TEMP_OUTDIR/cigar.txt $RAMA_OUTFILE
rm -r $RAMA_TEMP_OUTDIR

echo "aligning with minimap2"
MINIMAP_OUT=/private/groups/patenlab/mira/centrolign/analysis/variant_matrix/orthogonal_aligners/minimap2/${CHR}/
mkdir -p $MINIMAP_OUT

time docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
     mobinasri/long_read_aligner:v0.3.3 \
     minimap2 -x asm20 --eqx -c -f 0.05 --secondary=no -r 10000 \  
     -t 5 $FASTA1 $FASTA2 > $MINIMAP_OUT/${SAMP1}_${SAMP2}.paf
```
Convert minimap2 paf files to cigar strings
```sh
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")
for chr in "${chromosomes[@]}"
do

  ls /private/groups/patenlab/mira/centrolign/analysis/variant_matrix/orthogonal_aligners/minimap2/${chr}/ | while read line ; do
    s1=`echo ${line} | sed 's/\.paf//g' | cut -f1 -d"_"`
    s2=`echo ${line} | sed 's/\.paf//g' | cut -f2 -d"_"`

    # Extract CIGAR from primary alignment (tp:A:P), highest mapq
    grep 'tp:A:P' /private/groups/patenlab/mira/centrolign/analysis/variant_matrix/orthogonal_aligners/minimap2/${chr}/${line} | sort -k12,12nr | head -1 | grep -oP 'cg:Z:\K\S+' \
        > /private/groups/patenlab/mira/centrolign/analysis/variant_matrix/orthogonal_aligners/minimap2/${chr}/pairwise_cigar_${s1}_${s2}.txt
  done
done
```

Claude implemented Ti/Tv ratio calculator:
```sh
/private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/pairwise_titv_from_cigars.py

# unit tests
/private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/test_pairwise_titv_from_cigars.py
```

testing on chr5
```sh
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/pairwise_titv_from_cigars.py \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr5/subgroup_0/induced_pairwise_cigars/ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/chr5/ \
    chr5 \
    --indel_window 20 \
    --distance_csv /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/chr5_r2_QC_v2_centrolign_pairwise_distance.csv \
    --threshold 0.2 \
    > chr5_per_induced_pair_titv.tsv

#!/bin/sh
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/pairwise_titv_from_cigars.py \
    /private/groups/patenlab/mira/centrolign/analysis/variant_matrix/orthogonal_aligners/minimap2/chr5/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/chr5/ \
    chr5 \
    --indel_window 0 \
    --distance_csv /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/chr5_r2_QC_v2_centrolign_pairwise_distance.csv \
    --threshold 0.2 \
    > chr5_mm2_per_induced_pair_titv.tsv

#!/bin/sh
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/pairwise_titv_from_cigars.py \
    /private/groups/patenlab/mira/centrolign/analysis/variant_matrix/orthogonal_aligners/Rama/chr5/Rama_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/chr5/ \
    chr5 \
    --indel_window 0 \
    --distance_csv /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/chr5_r2_QC_v2_centrolign_pairwise_distance.csv \
    --threshold 0.2 \
    > chr5_rama_per_induced_pair_titv.tsv
```

Run on all chromosomes, all aligners: dist < 0.2, no indel based filtering
```sh
#!/bin/bash
#SBATCH --job-name=titv_orthogonal_aligners
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs_titv/%x.%j.log
#SBATCH --time=1:00:00
#SBATCH --array=[0]%24

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")
CHR=${chromosomes[$SLURM_ARRAY_TASK_ID]}

# Calculate for Centrolign
# get results for each subgroup
grep "release 2 QC v2" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/centrolign_results.csv | cut -f1,10 -d","  | grep -v "chr20," | grep "${CHR}," | while IFS=',' read -r subgroup cigar ; do
  python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/pairwise_titv_from_cigars.py \
      ${cigar}/pairwise_cigar_ \
      /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/${CHR}/ \
      ${CHR} \
      --indel_window 0 \
      --distance_csv /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/${CHR}_r2_QC_v2_centrolign_pairwise_distance.csv \
      --threshold 0.2 \
      &> /private/groups/patenlab/mira/centrolign/analysis/variant_matrix/orthogonal_aligners/pairwise_titv_from_cigars_results/dist_0.2/${subgroup}.centrolign.tsv
done

# Calculate for Unialigner
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/pairwise_titv_from_cigars.py \
    /private/groups/patenlab/mira/centrolign/analysis/variant_matrix/orthogonal_aligners/unialigner/${CHR}/unialigner_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/${CHR}/ \
    ${CHR} \
    --indel_window 0 \
    --distance_csv /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/${CHR}_r2_QC_v2_centrolign_pairwise_distance.csv \
    --threshold 0.2 \
    &> /private/groups/patenlab/mira/centrolign/analysis/variant_matrix/orthogonal_aligners/pairwise_titv_from_cigars_results/dist_0.2/${CHR}.unialigner.tsv

# Calculate for RAMA
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/pairwise_titv_from_cigars.py \
    /private/groups/patenlab/mira/centrolign/analysis/variant_matrix/orthogonal_aligners/Rama/${CHR}/Rama_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/${CHR}/ \
    ${CHR} \
    --indel_window 0 \
    --distance_csv /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/${CHR}_r2_QC_v2_centrolign_pairwise_distance.csv \
    --threshold 0.2 \
    &> /private/groups/patenlab/mira/centrolign/analysis/variant_matrix/orthogonal_aligners/pairwise_titv_from_cigars_results/dist_0.2/${CHR}.Rama.tsv

# Calculate for minimap2
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/pairwise_titv_from_cigars.py \
    /private/groups/patenlab/mira/centrolign/analysis/variant_matrix/orthogonal_aligners/minimap2/${CHR}/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/${CHR}/ \
    ${CHR} \
    --indel_window 0 \
    --distance_csv /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/${CHR}_r2_QC_v2_centrolign_pairwise_distance.csv \
    --threshold 0.2 \
    &> /private/groups/patenlab/mira/centrolign/analysis/variant_matrix/orthogonal_aligners/pairwise_titv_from_cigars_results/dist_0.2/${CHR}.minimap2.tsv
```


Run on all chromosomes, all aligners: Slurm script
```sh
#!/bin/bash
#SBATCH --job-name=titv_orthogonal_aligners
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs_titv/%x.%j.log
#SBATCH --time=1:00:00
#SBATCH --array=[0-23]%24

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")
CHR=${chromosomes[$SLURM_ARRAY_TASK_ID]}

DIST=$1
INDEL=$2
OUTDIR=/private/groups/patenlab/mira/centrolign/analysis/variant_matrix/orthogonal_aligners/pairwise_titv_from_cigars_results/dist_${DIST}_indel_${INDEL}

mkdir -p $OUTDIR
# Calculate for Centrolign induced
# get results for each subgroup
grep "release 2 QC v2" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/centrolign_results.csv | cut -f1,10 -d","  | grep -v "chr20," | grep "${CHR}" | while IFS=',' read -r subgroup cigar ; do
  python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/pairwise_titv_from_cigars.py \
      ${cigar}/pairwise_cigar_ \
      /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/${CHR}/ \
      ${CHR} \
      --indel_window ${INDEL} \
      --distance_csv /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/${CHR}_r2_QC_v2_centrolign_pairwise_distance.csv \
      --threshold ${DIST} \
      &> ${OUTDIR}/${subgroup}.centrolign_induced.tsv
done

# Calculate for Centrolign direct
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/pairwise_titv_from_cigars.py \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/${CHR}/pairwise_cigar/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/${CHR}/ \
    ${CHR} \
    --indel_window ${INDEL} \
    --distance_csv /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/${CHR}_r2_QC_v2_centrolign_pairwise_distance.csv \
    --threshold ${DIST} \
    &> ${OUTDIR}/${CHR}.centrolign_direct.tsv

# Calculate for Unialigner
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/pairwise_titv_from_cigars.py \
    /private/groups/patenlab/mira/centrolign/analysis/variant_matrix/orthogonal_aligners/unialigner/${CHR}/unialigner_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/${CHR}/ \
    ${CHR} \
    --indel_window ${INDEL} \
    --distance_csv /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/${CHR}_r2_QC_v2_centrolign_pairwise_distance.csv \
    --threshold ${DIST} \
    &> ${OUTDIR}/${CHR}.unialigner.tsv

# Calculate for RAMA
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/pairwise_titv_from_cigars.py \
    /private/groups/patenlab/mira/centrolign/analysis/variant_matrix/orthogonal_aligners/Rama/${CHR}/Rama_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/${CHR}/ \
    ${CHR} \
    --indel_window ${INDEL} \
    --distance_csv /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/${CHR}_r2_QC_v2_centrolign_pairwise_distance.csv \
    --threshold ${DIST} \
    &> ${OUTDIR}/${CHR}.Rama.tsv

# Calculate for minimap2
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/pairwise_titv_from_cigars.py \
    /private/groups/patenlab/mira/centrolign/analysis/variant_matrix/orthogonal_aligners/minimap2/${CHR}/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/${CHR}/ \
    ${CHR} \
    --indel_window ${INDEL} \
    --distance_csv /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/${CHR}_r2_QC_v2_centrolign_pairwise_distance.csv \
    --threshold ${DIST} \
    &> ${OUTDIR}/${CHR}.minimap2.tsv
```

Submitting slurm script:
```sh
cd /private/groups/patenlab/mira/centrolign/analysis/variant_matrix/orthogonal_aligners/pairwise_titv_from_cigars_results
mkdir -p logs_titv

# dist < 0.2, no indel filt
sbatch slurm.sh 0.2 0
# dist < 0.1, no indel filt
sbatch slurm.sh 0.1 0
# dist < 0.2, 20 bp indel filt
sbatch slurm.sh 0.2 20
# dist < 0.1, 20 bp indel filt
sbatch slurm.sh 0.1 20
```