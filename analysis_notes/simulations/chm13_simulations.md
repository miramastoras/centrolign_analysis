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
#!/bin/bash
#SBATCH --job-name=make-sim-cases
# Partition - This is the queue it goes in:
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=16gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-30
#SBATCH --output=logs/array_job_%A_task_%a.log
#SBATCH --time=6:00:00

date
hostname
pwd
echo "array job" $SLURM_ARRAY_TASK_ID

# note: sample size is handled in sbatch array size
CHR=12
DATE=20250324

SIMDIR=/private/groups/patenlab/mira/centrolign/simulations
OUTDIR=$SIMDIR/msa_chr"$CHR"_sim_cases_"$DATE"/case_"$SLURM_ARRAY_TASK_ID"/

GEN_TREE=/private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/generate_tree.py
SIM_CENTROMERE=/private/home/mmastora/progs/centrolign/build/sim_centromere

# the base array that we'll simulate
FASTA=/private/groups/patenlab/jeizenga/centromere/simulation/chm13_chr12_active_array.upper.fasta
BED=/private/groups/patenlab/jeizenga/centromere/simulation/chr12_shifted_hors.bed

# number of sequences in each case
N_SEQS=8
# expected number of generations to the root of the tree
N_GENS=200

mkdir -p $OUTDIR

# simulate the tree
source /private/groups/patenlab/jeizenga/centromere/venv/bin/activate
TREE=$OUTDIR/tree.txt
echo "making tree with" $N_SEQS "sequences and expected height" $N_GENS "generations"
time $GEN_TREE $N_SEQS $N_GENS > $TREE
# simulate the sequence
echo "generating sequences according to tree"
time $SIM_CENTROMERE -o $OUTDIR/sim -T $TREE $FASTA $BED
cat $OUTDIR/sim*.fasta > $OUTDIR/all_seqs.fasta
```

### 3. Run pairwise simulations and comparisons to other tools

Location:
```
/private/groups/patenlab/mira/centrolign/simulations/pairwise_simulations
```

Prepare simulated sequences
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
DATE=20250327

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
#SBATCH --partition=main
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=56gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=2-150
#SBATCH --output=logs/array_job_%A_task_%a.log
#SBATCH --time=1:00:00

date
hostname
pwd

# list cases in $SIMDIR/cases.txt
CHR=X
DATE=20231215

SIMDIR=/private/groups/patenlab/mira/centrolign/simulations/pairwise_simulations/
OUTPARDIR=$SIMDIR/pair_chr"$CHR"_sim_cases_"$DATE"/CASE=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$SIMDIR"/cases.txt)
CASEDIR=$SIMDIR/$CASE
WORKDIR=$SIMDIR/work

CENTROLIGN_OUTFILE=$CASEDIR/aln_centrolign.txt
UNIALIGNER_OUTFILE=$CASEDIR/aln_unialigner.txt

CENTROLIGN_DIR=/private/home/mmastora/progs/centrolign/build/
SCRIPTS_DIR=/private/groups/patenlab/mira/centrolign/github/centromere-scripts/

CENTROLIGN=$CENTROLIGN_DIR/centrolign
TRUTH_COMPARE=$CENTROLIGN_DIR/compare_truth_aln
UNIALIGNER=/private/groups/patenlab/jeizenga/GitHub/unialigner/tandem_aligner/build/bin/tandem_aligner
TO_RAW_SEQ=$SCRIPTS_DIR/data_processing_utils/fasta_to_raw_seq.py
ANALYZE_CASE=$SCRIPTS_DIR/benchmarking/analyze_pair_case.py

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
rm $TEMP_FASTA

echo "aligning with unaligner"
mkdir -p `dirname $UNIALIGNER_OUTFILE`
UNIALIGNER_TEMP_OUTDIR=$WORKDIR/tmp_out_"$SLURM_ARRAY_TASK_ID"
/usr/bin/time -v ${UNIALIGNER} --first $FASTA1 --second $FASTA2 -o $UNIALIGNER_TEMP_OUTDIR
mv $UNIALIGNER_TEMP_OUTDIR/cigar.txt $UNIALIGNER_OUTFILE
# delete the rest of the output
rm -r $UNIALIGNER_TEMP_OUTDIR

# do this from outside the directory to get more sensible output
cd $SIMDIR
$ANALYZE_CASE $CASE $TRUTH_COMPARE
```
Test out ramma
```
conda activate RaMA

RaMA -r sim_seq1.fasta -q sim_seq2.fasta -o /private/groups/patenlab/mira/centrolign/simulations/pairwise_simulations/pair_chr12_sim_cases_20250327/gen100/case_50/rama

/private/home/mmastora/progs/centrolign/build/centrolign -v 3 both.fasta > both.centrolign
```
```
RaMA -r /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/HG04184_mat/analysis/extract_hors_HPRC_outputs/HG04184_mat_HG04184.2_chr12_hor_array.fasta -q /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/HG00673_mat/analysis/extract_hors_HPRC_outputs/HG00673_mat_HG00673.2_chr12_hor_array.fasta -o /private/groups/patenlab/mira/centrolign/test/rama

RaMA -t 5 -r /private/groups/patenlab/mira/HG04184_mat_HG04184.2_chr12_hor_array.fasta -q /private/groups/patenlab/mira/HG00673_mat_HG00673.2_chr12_hor_array.fasta -o /private/groups/patenlab/mira/rama_test


cat /private/groups/patenlab/mira/HG04184_mat_HG04184.2_chr12_hor_array.fasta /private/groups/patenlab/mira/HG00673_mat_HG00673.2_chr12_hor_array.fasta > /private/groups/patenlab/mira/centrolign/test/hprc_test.fa
/private/home/mmastora/progs/centrolign/build/centrolign -v 3 /private/groups/patenlab/mira/centrolign/test/hprc_test.fa > /private/groups/patenlab/mira/centrolign/test/hprc_test_centrolign

```
