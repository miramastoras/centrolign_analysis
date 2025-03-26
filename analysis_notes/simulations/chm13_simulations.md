## Simulating HOR sequences to test and evaluate centrolign


### 1. Identify cyclic HOR arrays without SVs

https://pmc.ncbi.nlm.nih.gov/articles/PMC9248890/

```
chr 2,3,4,6,7,10,11,12,14,15,16,17,20,21,22,X,Y
```

AS_HOR bed from Hailey for CHM13: https://public.gi.ucsc.edu/~hloucks/CenSat/CHM13/chm13v2.0.labels.as_hor.bed

Generate CHM13 files for all of these chromosomes
```

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
DATE=20231215

CENTRO_DIR=/private/groups/patenlab/jeizenga/centromere
SIMDIR=$CENTRO_DIR/simulation/
OUTPARDIR=$SIMDIR/pair_chr"$CHR"_sim_cases_"$DATE"/

SIM_CENTROMERE=/private/groups/patenlab/jeizenga/GitHub/centrolign/build/sim_centromere-429efbd

# the base array that we'll simulate
FASTA=/private/groups/patenlab/jeizenga/centromere/simulation/chm13_chr12_active_array.upper.fasta
BED=/private/groups/patenlab/jeizenga/centromere/simulation/chr12_shifted_hors.bed

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
