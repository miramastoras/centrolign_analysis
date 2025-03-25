## Simulating HOR sequences to test and evaluate centrolign

### 1. Identify cyclic HOR arrays without SVs

https://pmc.ncbi.nlm.nih.gov/articles/PMC9248890/

chr 2,3,4,6,7,10,11,12,14,15,16,17,20,21,22,X,Y

### 2. make sim cases

```
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
