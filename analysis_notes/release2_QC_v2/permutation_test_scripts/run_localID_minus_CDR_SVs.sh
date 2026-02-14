#!/bin/bash
#SBATCH --job-name=localID_minus_CDR_SVs
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=logs/%x.%j.log
#SBATCH --time=7-00:00

source /private/home/mmastora/miniconda3/etc/profile.d/conda.sh
conda activate ipynb 

python /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/permutation_test_scripts/localID_minus_CDR_permutation_per_chrom_SVs.py \
    --n-reps 10000 \
    --n-workers 24 \
    --output-dir /private/groups/patenlab/mira/centrolign/analysis/local_identity/permutation/localID_minus_CDR_SVs_10k
