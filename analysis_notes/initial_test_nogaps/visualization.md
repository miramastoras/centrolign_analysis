### Chr 12 pairwise heatmaps

Run centrolign outputting pairwise cigars
```
#!/bin/bash
#SBATCH --job-name=centrolign_chr12_nogaps_all_samples_pairwise
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=700gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=centrolign_restart2_12_12_24.log
#SBATCH --time=7-00:00

/private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/jobstore/ \
    -T /private/groups/patenlab/jeizenga/centromere/chr12/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.nwk.txt \
    -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/pairwise_cigars/ \
    -R \
    --threads 32 \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/initial_test_no_gaps_chr12.fasta \
    > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/initial_test_no_gaps_chr12.centrolign.gfa
```

Convert pairwise cigars to distances for tree heatmap
```
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12

python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/cigar_to_distance.py /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/pairwise_cigars/
```

Plot tree heatmap
```
python3 pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/tree_heatmap/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.nwk.txt \
        -s /Users/miramastoras/Desktop/tree_heatmap_chr12/samples.txt \
        -p /Users/miramastoras/Desktop/tree_heatmap_chr12/pairwise_distance.csv \
        -o /Users/miramastoras/Desktop/tree_heatmap_chr12/
```
Plot just trio samples
```
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/tree_heatmap/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.nwk.txt \
        -s /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.txt \
        -p /Users/miramastoras/Desktop/tree_heatmap_chr12/pairwise_distance.csv \
        -o /Users/miramastoras/Desktop/plot_all_tree_methods_just_trios/cenhap_tree_centrolign_MSA_distances_chr12_trios_
```
