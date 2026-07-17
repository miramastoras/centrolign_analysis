### Re-running chrY for the censat paper with new parameters 

1. Pull and rebuild centrolign 

2. Submit all pairs with config files 

Run centrolign all pairs - Chr Y
```sh
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chrY_072326/logs

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chrY_072326

sbatch \
    --job-name=chrY_r2QCv2_all_pairs_rerun \
    --array=[1-990]%128 \
    --export=CHR=chrY \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/all_pairs_chrY.sh
```

3. Create distance matrices 

```sh
GITHUB=/private/groups/patenlab/mira/centrolign/github/centrolign_analysis/
BATCH_DIR=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2

# Get distance matrix from all pairs
python3 ${GITHUB}/scripts/cigar_to_distance.py -d 1 \
-a ${BATCH_DIR}/all_pairs/chrY_072326/pairwise_cigar/ \
-o ${BATCH_DIR}/all_pairs/distance_matrices/chrY_072326_r2_QC_v2_centrolign_

# infer neighbor-joining tree from distance matrix
docker run -u `id -u`:`id -g` -v /private/groups:/private/groups/ \
    miramastoras/centromere_scripts:v0.1.4 \
    python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/data_exploration/infer_tree.py \
    ${BATCH_DIR}/all_pairs/chrY_072326/pairwise_cigar/ \
    > ${BATCH_DIR}/all_pairs/nj_trees/chrY_072326_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk

# convert tree to format needed for plotting
docker run -u `id -u`:`id -g` -v /private/groups:/private/groups/ \
    miramastoras/centromere_scripts:v0.1.4 \
    python3 ${GITHUB}/scripts/convert_tree_format.py \
    -t ${BATCH_DIR}/all_pairs/nj_trees/chrY_072326_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk

```
4. Submit MSA

Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrY_072326/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrY_072326/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrY_072326/
```


```sh
#!/bin/bash
#SBATCH --job-name=chrY_072326_MSA
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=800gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=7-00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  --config /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrY_072326/chrY_config.yaml \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrY_072326/chrY.centrolign.gfa
```

5. Submit induced pairwise 

```sh
#!/bin/bash
#SBATCH --job-name=chrY_072326_MSA
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=800gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=7-00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  --config /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrY_072326/chrY_induced_pairwise_config.yaml \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrY_072326/chrY.centrolign.gfa
```

6. Get distance matrix for induced pairwise cigar
```sh
GITHUB=/private/groups/patenlab/mira/centrolign/github/centrolign_analysis/
# Get distance matrix from all pairs
python3 ${GITHUB}/scripts/cigar_to_distance.py -d 1 \
-a /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrY_072326/induced_pairwise_cigars/ \
-o /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrY_072326/chrY_072326_r2_QC_v2_centrolign_induced_
```

7. Plot heatmap
```sh
chromosomes=("chrY")

for chr in "${chromosomes[@]}"
do
  python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap_v2.py \
    -t /Users/miramastoras/Desktop/chrY_072326_r2_QC_v2_centrolign_all_pairs_nj_tree.format5.nwk \
    -s /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}.samples.txt \
    -p /Users/miramastoras/Desktop/chrY_072326_r2_QC_v2_centrolign_pairwise_distance.csv \
    -m "Centrolign all pairs distances" \
    -n "${chr} NJ tree" \
    -d "All pairs Distances" \
    -o /Users/miramastoras/Desktop/chrY_072326_ --no_labels
done
```