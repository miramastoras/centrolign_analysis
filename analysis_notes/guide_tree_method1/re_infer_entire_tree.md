## Testing reinferring the entire guide tree using centrolign pairwise distances and neighbor-joining

### 1. Run infer_tree on pairwise cigar strings for all samples

```
# infer tree requires pairwise text files have aln or cigar in their title
ls | while read line ; do mv ${line} pairwise_cigar${line} ; done

conda activate skbio

python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/data_exploration/infer_tree.py /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/pairwise_cigars/ > /private/groups/patenlab/mira/centrolign/guide_tree_testing/method1/re_infer_entire_tree/chr12_inferred_tree.nwk
```
Convert 0 branch lengths and scale up
```
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/scale_branch_lengths.py \
    -t /Users/miramastoras/Desktop/tree_heatmap/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.nwk.txt \
    -n /Users/miramastoras/Desktop/re_infer_whole_tree/chr12_inferred_tree.nwk \
    -i /Users/miramastoras/Desktop/tree_heatmap_chr12/samples.txt \
    -o /Users/miramastoras/Desktop/re_infer_whole_tree/chr12_inferred_tree.rescaled.nwk.txt
```

### 2. Rerun pairwise heatmap script and tanglegram

```
conda activate tree_python

python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/re_infer_whole_tree/chr12_inferred_tree.rescaled.nwk \
        -s /Users/miramastoras/Desktop/tree_heatmap_chr12/samples.txt  \
        -p /Users/miramastoras/Desktop/tree_heatmap_chr12/pairwise_distance.csv \
        -o /Users/miramastoras/Desktop/re_infer_whole_tree/re_infer_whole_tree
```

### 3. Restart centrolign

First submit centrolign in MSA mode. Then restart from last subproblem outputting pairwise alignments in multithreaded mode
```
#!/bin/bash
#SBATCH --job-name=centrolign_chr12_new_tree
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=700gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign.log
#SBATCH --time=7-00:00

/private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S /private/groups/patenlab/mira/centrolign/guide_tree_testing/method1/re_infer_entire_tree/jobstore/ \
    -T /private/groups/patenlab/mira/centrolign/guide_tree_testing/method1/re_infer_entire_tree/chr12_inferred_tree.rescaled.nwk.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/initial_test_no_gaps_chr12.fasta \
    > /private/groups/patenlab/mira/centrolign/guide_tree_testing/method1/re_infer_entire_tree/initial_test_no_gaps_chr12.re_infer_whole_tree.centrolign.gfa
```
