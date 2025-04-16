## Simulating HOR sequences to test and evaluate centrolign


### 1. Identify cyclic HOR arrays without SVs

https://pmc.ncbi.nlm.nih.gov/articles/PMC9248890/

```
chr 2,3,4,6,7,10,11,12,14,15,16,17,20,21,22,X,Y
```

AS_HOR bed from Hailey for CHM13: https://public.gi.ucsc.edu/~hloucks/CenSat/CHM13/chm13v2.0.labels.as_hor.bed

Extract HOR from chm13 for each chromosome
```sh
chromosomes=("chr2" "chr3" "chr4" "chr6" "chr7" "chr10" "chr11" "chr12" "chr14" "chr15" "chr16" "chr17" "chr20" "chr21" "chr22" "chrX" "chrY")

mkdir -p per_chrom/work

for chr in "${chromosomes[@]}"
do
  # bedtools intersect alpha annotation with censat annotation to just get active hors
  echo "processing ${chr} "
  grep -w $chr chm13v2.0_censat_v2.1.bed | grep "hor" | grep -v "dhor" > per_chrom/work/chm13v2.0_censat_v2.1.${chr}.hor.bed

  bedtools intersect -wa -a chm13v2.0.labels.as_hor.bed -b per_chrom/work/chm13v2.0_censat_v2.1.${chr}.hor.bed > per_chrom/work/chm13v2.0.labels.as_hor.${chr}.active.bed

  # get start and end of active array from as_hor bed
  awk 'NR==1{start=$2; chrom=$1} END{print chrom"\t"start"\t"$3}' per_chrom/work/chm13v2.0.labels.as_hor.${chr}.active.bed > per_chrom/work/chm13v2.0.labels.as_hor.${chr}.active.mrg.bed

  # extract fasta sequence for active hor array
  bedtools getfasta -fi /private/groups/patenlab/mira/data/chm13v2.0.fa -bed per_chrom/work/chm13v2.0.labels.as_hor.${chr}.active.mrg.bed | sed '/^>/b; s/[a-z]/\U&/g' | sed '/^>/s/:.*//' > per_chrom/chm13v2.0.${chr}.active_hor.upper.fa

  # shift hor arrays over to start at 0
  awk 'NR==1 {shift=$2} {print $1"\t"($2-shift)"\t"($3-shift)"\t"$4}' per_chrom/work/chm13v2.0.labels.as_hor.${chr}.active.bed > per_chrom/chm13v2.0.labels.as_hor.${chr}.active.shifted.bed

done
```

### 2. Run MSA simulations


#### Step 1: make sim cases

```sh

mkdir -p /private/groups/patenlab/mira/centrolign/simulations/MSA_simulations/make_sim_cases_slurm_logs

git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/ pull

cd /private/groups/patenlab/mira/centrolign/simulations/MSA_simulations

chromosomes=("chr2" "chr3" "chr4" "chr6" "chr7" "chr10" "chr11" "chr12" "chr14" "chr15" "chr16" "chr17" "chr20" "chr21" "chr22" "chrX" "chrY")

for chr in "${chromosomes[@]}"
do
    sbatch /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/simulations/slurm_scripts/make_sim_cases_MSA.sh $chr
done
```

#### Step 2: Run centrolign MSA on simulated seqs and analyze results

```sh

mkdir -p /private/groups/patenlab/mira/centrolign/simulations/MSA_simulations/simulated_centrolign_slurm_logs

git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/ pull

cd /private/groups/patenlab/mira/centrolign/simulations/MSA_simulations

chromosomes=("chr2" "chr3" "chr4" "chr6" "chr7" "chr10" "chr11" "chr12" "chr14" "chr15" "chr16" "chr17" "chr20" "chr21" "chr22" "chrX" "chrY")

for chr in "${chromosomes[@]}"
do
    sbatch /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/simulations/slurm_scripts/simulated_centrolign_MSA.sh $chr
done
```
#### Plot results for MSA simulations

```sh
cd /private/groups/patenlab/mira/centrolign/simulations/MSA_simulations

chromosomes=("chr2" "chr3" "chr4" "chr6" "chr7" "chr10" "chr11" "chr12" "chr14" "chr15" "chr16" "chr17" "chr20" "chr21" "chr22" "chrX" "chrY")

mkdir -p /private/groups/patenlab/mira/centrolign/simulations/MSA_simulations/summary_tables

for chr in "${chromosomes[@]}"
do
    echo "Processing $chr"
    cat msa_${chr}_sim_cases_20250402/*/aln_summary_table.txt > summary_tables/msa_${chr}_sim_cases_20250402_aln_summary_tables.txt
done
```
Plot results
```sh
cd /private/groups/patenlab/mira/centrolign/simulations/MSA_simulations/summary_tables

ls | while read line ; do
    chr=`echo $line | cut -f2 -d"_"`
    Rscript /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/msa_simulations.R $line $chr /private/groups/patenlab/mira/centrolign/simulations/MSA_simulations/png_plots/${chr}_MSA_simulations
  done
```

## 3. Run pairwise simulations and comparisons to other tools

Location:
```
/private/groups/patenlab/mira/centrolign/simulations/pairwise_simulations
```

#### Step 1: make sim cases

```sh

mkdir -p /private/groups/patenlab/mira/centrolign/simulations/pairwise_simulations/make_sim_cases_slurm_logs

git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/ pull

cd /private/groups/patenlab/mira/centrolign/simulations/pairwise_simulations

chromosomes=("chr2" "chr3" "chr4" "chr6" "chr7" "chr10" "chr11" "chr12" "chr14" "chr15" "chr16" "chr17" "chr20" "chr21" "chr22" "chrX" "chrY")

for chr in "${chromosomes[@]}"
do
    sbatch /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/simulations/slurm_scripts/make_sim_cases_pairwise.sh $chr
done
```
#### Step 2: run centrolign and other pairwise aligners
```sh

mkdir -p /private/groups/patenlab/mira/centrolign/simulations/pairwise_simulations/simulated_centrolign_slurm_logs

chromosomes=("chr2" "chr3" "chr4" "chr6" "chr7" "chr10" "chr11" "chr12" "chr14" "chr15" "chr16" "chr17" "chr20" "chr21" "chr22" "chrX" "chrY")

git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/ pull

cd /private/groups/patenlab/mira/centrolign/simulations/pairwise_simulations

for chr in "${chromosomes[@]}"
do
    sbatch /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/simulations/slurm_scripts/simulated_centrolign_pairwise.sh $chr
done
```


```sh
#!/bin/bash
# Slurm script to run the sim_centromere script to generate pairwise sequence
# alignment problems
#SBATCH --job-name=make-pair-sim-cases
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=3gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-60
#SBATCH --output=logs/array_job_%A_task_%a.log
#SBATCH --time=1:00:00

date
hostname
pwd
echo "array job" $SLURM_ARRAY_TASK_ID

# note: sample size is handled in sbatch array size
CHR=12
DATE=20250331

SIMDIR=/private/groups/patenlab/mira/centrolign/simulations/pairwise_simulations/
OUTPARDIR=$SIMDIR/pair_chr"$CHR"_sim_cases_"$DATE"/

SIM_CENTROMERE=/private/home/mmastora/progs/centrolign/build/sim_centromere

# the base array that we'll simulate
FASTA=/private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom/chm13v2.0.chr${CHR}.active_hor.upper.fa
BED=/private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom/chm13v2.0.labels.as_hor.chr${CHR}.active.shifted.bed

GENS=("25" "50" "100" "150" "200" "300")
i=$(($SLURM_ARRAY_TASK_ID % 6))
GEN=${GENS[$i]}
OUTDIR=$OUTPARDIR/gen"$GEN"/case_"$SLURM_ARRAY_TASK_ID"/
mkdir -p $OUTDIR
# simulate the sequence
echo "generating sequence pair" $SLURM_ARRAY_TASK_ID "with" $GEN "generations"
/usr/bin/time -v $SIM_CENTROMERE -g $GEN -o $OUTDIR/sim $FASTA $BED
# record this case for the centrolign file
echo gen"$GEN"/case_"$SLURM_ARRAY_TASK_ID" >> $OUTPARDIR/cases.txt

# this is silly, but it makes the evaluation logic much simpler if we always have a tree
printf "(seq1,seq2):%d;\n" $GEN > $OUTDIR/tree.txt
cat $OUTDIR/sim*.fasta > $OUTDIR/all_seqs.fasta
```

Run centrolign and the alternative tools
```sh
#!/bin/bash
# Slurm script to align simulated sequences pairwise with centrolign, Unialigner, winnowmap, and RaMa
# and to analyze the accuracy, storing the results in a table
#SBATCH --job-name=simulated-centrolign
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=56gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1
#SBATCH --output=logs/analysis_array_job_%A_task_%a.log
#SBATCH --time=12:00:00

date
hostname
pwd

# list cases in $SIMDIR/cases.txt
CHR=12
DATE=20250331

SIMDIR=/private/groups/patenlab/mira/centrolign/simulations/pairwise_simulations/pair_chr"$CHR"_sim_cases_"$DATE"/
CASE=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$SIMDIR"/cases.txt)
CASEDIR=$SIMDIR/$CASE
WORKDIR=$SIMDIR/work

CENTROLIGN_OUTFILE=$CASEDIR/aln_centrolign.txt
UNIALIGNER_OUTFILE=$CASEDIR/aln_unialigner.txt
RAMA_OUTFILE=$CASEDIR/aln_rama.txt
WFA_OUTFILE=$CASEDIR/aln_wfa.txt
#WINNOWMAP_OUTFILE=$CASEDIR/aln_winnowmap.txt
#MINIMAP2_OUTFILE=$CASEDIR/aln_minimap2.txt

WFA=/private/groups/patenlab/jeizenga/GitHub/WFA2-lib/bin/align_benchmark
CENTROLIGN_DIR=/private/home/mmastora/progs/centrolign/build/
SCRIPTS_DIR=/private/groups/patenlab/mira/centrolign/github/centromere-scripts/

# scores: M,X,O1,E1,O2,E2
# heuristic: min distance, max distance from lead, reduction interval
WFA_PARAMS="-a gap-affine2p-wfa --affine2p-penalties -20,80,100,30,5000,1 --wfa-heuristic wfa-adaptive --wfa-heuristic-parameters 1000,20000,50 --wfa-memory ultralow --wfa-span global"


CENTROLIGN=$CENTROLIGN_DIR/centrolign
TRUTH_COMPARE=$CENTROLIGN_DIR/compare_truth_aln
UNIALIGNER=/private/groups/patenlab/mira/centrolign/github/unialigner/tandem_aligner/build/bin/tandem_aligner
TO_RAW_SEQ=$SCRIPTS_DIR/data_processing_utils/fasta_to_raw_seq.py
ANALYZE_CASE=/private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/analyze_pair_case.py

mkdir -p $WORKDIR
cd $WORKDIR

echo "simulation:" $CASEDIR

FASTA1=${CASEDIR}/sim_seq1.fasta
FASTA2=${CASEDIR}/sim_seq2.fasta

echo "fasta 1:" $FASTA1
echo "fasta 2:" $FASTA2

echo "aligning with centrolign"
TEMP_FASTA=${WORKDIR}/sim_joined_"$SLURM_ARRAY_TASK_ID".fa
mkdir -p `dirname $TEMP_FASTA`
cat $FASTA1 $FASTA2 > $TEMP_FASTA
mkdir -p `dirname $CENTROLIGN_OUTFILE`
/usr/bin/time -v ${CENTROLIGN} -v 3 $TEMP_FASTA > $CENTROLIGN_OUTFILE
#rm $TEMP_FASTA

echo "aligning with unaligner"
mkdir -p `dirname $UNIALIGNER_OUTFILE`
UNIALIGNER_TEMP_OUTDIR=$WORKDIR/uni_tmp_out_"$SLURM_ARRAY_TASK_ID"
time ${UNIALIGNER} --first $FASTA1 --second $FASTA2 -o $UNIALIGNER_TEMP_OUTDIR
mv $UNIALIGNER_TEMP_OUTDIR/cigar.txt $UNIALIGNER_OUTFILE
# delete the rest of the output
#rm -r $UNIALIGNER_TEMP_OUTDIR

echo "aligning with RAMA"
RAMA_TEMP_OUTDIR=$WORKDIR/rama_tmp_out_"$SLURM_ARRAY_TASK_ID"
mkdir -p `dirname $RAMA_OUTFILE`
time docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
    miramastoras/rama:latest ./RaMA -t 5 \
    -r $FASTA1 \
    -q $FASTA2 \
    -o $RAMA_TEMP_OUTDIR

mv $RAMA_TEMP_OUTDIR/cigar.txt $RAMA_OUTFILE

# delete the rest of the output
rm -r $RAMA_TEMP_OUTDIR

echo "aligning with WFA"
RAW_SEQ_TEMP=${WORKDIR}/tmp_raw_seq_"$SLURM_ARRAY_TASK_ID".txt
WFA_TEMP_OUT=${WORKDIR}/tmp_wfa_out_"$SLURM_ARRAY_TASK_ID".txt
$TO_RAW_SEQ $FASTA1 > $RAW_SEQ_TEMP
$TO_RAW_SEQ $FASTA2 >> $RAW_SEQ_TEMP
# limit memory to 32 gB and runtime to 30 min, but don't consider it a failure if we don't get it
true || timeout -v 30m ulimit -m 33554432 $WFA $WFA_PARAMS -i $RAW_SEQ_TEMP -o $WFA_TEMP_OUT
# remove the score from the output
if [ -f $WFA_TEMP_OUT ]; then
   cut -f 2 $WFA_TEMP_OUT > $WFA_OUTFILE
else
   touch $WFA_OUTFILE
fi
rm -f $RAW_SEQ_TEMP
rm -f $WFA_TEMP_OUT

# echo "aligning with winnowmap"
#
# WINNOWMAP_TEMP_OUTDIR=$WORKDIR/winnow_tmp_out_"$SLURM_ARRAY_TASK_ID"
#
# mkdir -p ${WINNOWMAP_TEMP_OUTDIR}
# mkdir -p ${WINNOWMAP_OUTFILE}
#
# time docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
#     mobinasri/long_read_aligner:v0.3.3 \
#     meryl count k=19 output ${WINNOWMAP_TEMP_OUTDIR}/merylDB $FASTA1
#
# time docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
#     mobinasri/long_read_aligner:v0.3.3 \
#     meryl print greater-than distinct=0.9998 \
#     ${WINNOWMAP_TEMP_OUTDIR}/merylDB \
#     > ${WINNOWMAP_TEMP_OUTDIR}/repetitive_k19.txt
#
# docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
#     mobinasri/long_read_aligner:v0.3.3 \
#     winnowmap \
#     -W ${WINNOWMAP_TEMP_OUTDIR}/repetitive_k19.txt \
#     --eqx -cx asm20 -t 1 \
#     $FASTA1 $FASTA2 \
#     > ${WINNOWMAP_OUTFILE}
#
# #rm -r $WINNOWMAP_TEMP_OUTDIR

# do this from outside the directory to get more sensible output
cd $SIMDIR
python3 $ANALYZE_CASE $CASE $TRUTH_COMPARE
```
