### Using centrolign all pairs to create HOR trees as inputs

Generate NJ trees for all samples (except chr6 and 12 since those have been created)
```sh
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr7" "chr8" "chr9" "chr10" "chr11" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

for chr in "${chromosomes[@]}"
do
    docker run -u `id -u`:`id -g` -v /private/groups:/private/groups/ \
        miramastoras/centromere_scripts:v0.1.4 \
        python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/data_exploration/infer_tree.py \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/${chr}/pairwise_cigar/ \
        > /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_all_pairs/${chr}_r2_centrolign_all_pairs_nj_tree.nwk

    docker run -u `id -u`:`id -g` -v /private/groups:/private/groups/ \
        miramastoras/centromere_scripts:v0.1.4 \
        python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/convert_tree_format.py \
        -t /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_all_pairs/${chr}_r2_centrolign_all_pairs_nj_tree.nwk
  done
```
```sh
docker run -u `id -u`:`id -g` -v /private/groups:/private/groups/ \
    miramastoras/centromere_scripts:v0.1.4 \
    python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/convert_tree_format.py \
    -t /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_all_pairs/test.nwk \
    -a 1 -b 5
```

#### April 28, 2025 - Submit chromosomes with lowest sample size

chrY, chr3, chr4, chr18, chr17  
```
398 chr1
356 chr10
422 chr11
385 chr12
303 chr13
348 chr14
331 chr15
402 chr16
139 chr17
130 chr18
392 chr19
433 chr2
343 chr20
275 chr21
333 chr22
160 chr3
196 chr4
375 chr5
191 chr6
337 chr7
358 chr8
407 chr9
315 chrX
 49 chrY
```
Run centrolign on chrY
```sh
#!/bin/bash
#SBATCH --job-name=centrolign_release2_chrY_allpairs_NJ_HOR_tree
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=700gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign_%x.%j.log
#SBATCH --time=7-00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/MSA_HOR_all_pairs_NJ_tree/chrY/jobstore/ \
    -T /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_all_pairs/chrY_r2_centrolign_all_pairs_nj_tree.format5.nwk \
    /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/centrolign_fastas/HPRC_release2_HORs_chrY.fasta \
    > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/MSA_HOR_all_pairs_NJ_tree/chrY/HPRC_r2.chrY.allpairs_HOR_NJ.centrolign.gfa
```
Run centrolign on chr17
```sh
#!/bin/bash
#SBATCH --job-name=centrolign_release2_chr17_allpairs_NJ_HOR_tree
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=700gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign_%x.%j.log
#SBATCH --time=7-00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/MSA_HOR_all_pairs_NJ_tree/chr17/jobstore/ \
    -T /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_all_pairs/chr17_r2_centrolign_all_pairs_nj_tree.format5.nwk \
    /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/centrolign_fastas/HPRC_release2_HORs_chr17.fasta \
    > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/MSA_HOR_all_pairs_NJ_tree/chr17/HPRC_r2.chr17.allpairs_HOR_NJ.centrolign.gfa
```

Run centrolign on chr4
```sh
#!/bin/bash
#SBATCH --job-name=centrolign_release2_chr4_allpairs_NJ_HOR_tree
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=700gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign_%x.%j.log
#SBATCH --time=7-00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/MSA_HOR_all_pairs_NJ_tree/chr4/jobstore/ \
    -T /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_all_pairs/chr4_r2_centrolign_all_pairs_nj_tree.format5.nwk \
    /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/centrolign_fastas/HPRC_release2_HORs_chr4.fasta \
    > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/MSA_HOR_all_pairs_NJ_tree/chr4/HPRC_r2.chr4.allpairs_HOR_NJ.centrolign.gfa
```
