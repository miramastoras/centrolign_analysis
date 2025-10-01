#!/bin/bash
# Slurm script to use centrolign as a pairwise aligner on all pairs of sequences
# from an input directory

#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=56gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/array_job_%A_task_%a.log
#SBATCH --time=12:00:00

GITHUB=/private/groups/patenlab/mira/centrolign/github/centrolign_analysis/
BATCH_DIR=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2

# Get distance matrix from all pairs
python3 ${GITHUB}/scripts/cigar_to_distance.py -d 1 \
-a ${BATCH_DIR}/all_pairs/${CHR}/pairwise_cigar/ \
-o ${BATCH_DIR}/all_pairs/distance_matrices/${CHR}_r2_QC_v2_centrolign_

# infer neighbor-joining tree from distance matrix
docker run -u `id -u`:`id -g` -v /private/groups:/private/groups/ \
    miramastoras/centromere_scripts:v0.1.4 \
    python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/data_exploration/infer_tree.py \
    ${BATCH_DIR}/all_pairs/${CHR}/pairwise_cigar/ \
    > ${BATCH_DIR}/all_pairs/nj_trees/${CHR}_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk

# convert tree to format needed for plotting
docker run -u `id -u`:`id -g` -v /private/groups:/private/groups/ \
    miramastoras/centromere_scripts:v0.1.4 \
    python3 ${GITHUB}/scripts/convert_tree_format.py \
    -t ${BATCH_DIR}/all_pairs/nj_trees/${CHR}_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk
