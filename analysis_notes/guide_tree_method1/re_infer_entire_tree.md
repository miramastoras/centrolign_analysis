## Testing reinferring the entire guide tree using centrolign pairwise distances and neighbor-joining

### 1. Run infer_tree on pairwise cigar strings for all samples

```
# infer tree requires pairwise text files have aln or cigar in their title
ls | while read line ; do mv ${line} pairwise_cigar${line} ; done

conda activate skbio

python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/data_exploration/infer_tree.py /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/pairwise_cigars/ > /private/groups/patenlab/mira/centrolign/guide_tree_testing/method1/re_infer_entire_tree/chr12_inferred_tree.nwk
```
Convert 0 branch lengths and scale up 

### 2. Rerun pairwise heatmap script

```
conda activate tree_python

python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/re_infer_whole_tree/chr12_inferred_tree.nwk \
        -s /Users/miramastoras/Desktop/tree_heatmap_chr12/samples.txt  \
        -p /Users/miramastoras/Desktop/tree_heatmap_chr12/pairwise_distance.csv \
        -o /Users/miramastoras/Desktop/re_infer_whole_tree/re_infer_whole_tree
```

### 3. Restart centrolign
