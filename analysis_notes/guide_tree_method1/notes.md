### Testing new guide tree method - chr12 initial test no gaps

In this method we freeze the deep nodes of the tree, then re-infer the subtrees using neighbor joining and the pairwise distances from centrolign. Starting off we'll test splitting the tree into 5 subgroups, and using chr12.

#### 1. Plot pairwise alignments for centrolign on 128 HPRC samples, chr12

![chr12_tree](pics/chr12_initial_test_nogaps_pairwise_tree_heatmap.png)

```
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/tree_heatmap/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.nwk.txt \
        -s /Users/miramastoras/Desktop/tree_heatmap_chr12/samples.txt \
        -p /Users/miramastoras/Desktop/tree_heatmap_chr12/pairwise_distance.csv \
        -o /Users/miramastoras/Desktop/tree_heatmap_chr12/
```

#### 2. Split fasta and tree into 5 groups using [split_fasta_by_tree](https://github.com/jeizenga/centromere-scripts/blob/main/benchmarking/split_fasta_by_tree.py)

```
conda create -n skbio
conda install anaconda::scikit-bio
conda activate skbio

python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/split_fasta_by_tree.py \
    -t /private/groups/patenlab/jeizenga/centromere/chr12/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.nwk.txt \
    -f /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/initial_test_no_gaps_chr12.fasta \
    -n 5 \
    -o /private/groups/patenlab/mira/centrolign/guide_tree_testing/method1/chr12_initial_test_nogaps/split_fasta_by_tree/
```
#### 3. Restart from subproblem for each of the 5 groups, outputting pairwise alignments

```
#!/bin/bash
#SBATCH --job-name=centrolign_subgroup_4_chr12
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=700gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=centrolign_subgroup_4.log
#SBATCH --time=1:00:00
#SBATCH --array=0-5

/private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/jobstore/ \
    -T /private/groups/patenlab/mira/centrolign/guide_tree_testing/method1/chr12_initial_test_nogaps/split_fasta_by_tree/subgroup_4_tree.nwk \
    -A /private/groups/patenlab/mira/centrolign/guide_tree_testing/method1/chr12_initial_test_nogaps/centrolign_restart_subgroups/subgroup_4/pairwise_cigars/ \
    -R \
    --threads 32 \
    /private/groups/patenlab/mira/centrolign/guide_tree_testing/method1/chr12_initial_test_nogaps/split_fasta_by_tree/subgroup_4_seqs.fasta \
    > /private/groups/patenlab/mira/centrolign/guide_tree_testing/method1/chr12_initial_test_nogaps/centrolign_restart_subgroups/subgroup_4/subgroup_4_chr12.centrolign.gfa
```
subgroup 3 has just one sample in it, HG04115, so skipping that in future steps

#### 4. infer tree for each of the subgroups

```
# infer tree requires pairwise text files have aln or cigar in their title
ls | while read line ; do mv ${line} pairwise_cigar${line} ; done

conda activate skbio

python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/data_exploration/infer_tree.py /private/groups/patenlab/mira/centrolign/guide_tree_testing/method1/chr12_initial_test_nogaps/centrolign_restart_subgroups/subgroup_2/pairwise_cigars/ > /private/groups/patenlab/mira/centrolign/guide_tree_testing/method1/chr12_initial_test_nogaps/centrolign_restart_subgroups/subgroup_2/infer_tree/subgroup2_inferred_tree.nwk

python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/data_exploration/infer_tree.py /private/groups/patenlab/mira/centrolign/guide_tree_testing/method1/chr12_initial_test_nogaps/centrolign_restart_subgroups/subgroup_4/pairwise_cigars/ > /private/groups/patenlab/mira/centrolign/guide_tree_testing/method1/chr12_initial_test_nogaps/centrolign_restart_subgroups/subgroup_4/infer_tree/subgroup4_inferred_tree.nwk

python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/data_exploration/infer_tree.py /private/groups/patenlab/mira/centrolign/guide_tree_testing/method1/chr12_initial_test_nogaps/centrolign_restart_subgroups/subgroup_0/pairwise_cigars/ > /private/groups/patenlab/mira/centrolign/guide_tree_testing/method1/chr12_initial_test_nogaps/centrolign_restart_subgroups/subgroup_0/infer_tree/subgroup0_inferred_tree.nwk

python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/data_exploration/infer_tree.py /private/groups/patenlab/mira/centrolign/guide_tree_testing/method1/chr12_initial_test_nogaps/centrolign_restart_subgroups/subgroup_2/pairwise_cigars/ > /private/groups/patenlab/mira/centrolign/guide_tree_testing/method1/chr12_initial_test_nogaps/centrolign_restart_subgroups/subgroup_2/infer_tree/subgroup2_inferred_tree.nwk

python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/data_exploration/infer_tree.py /private/groups/patenlab/mira/centrolign/guide_tree_testing/method1/chr12_initial_test_nogaps/centrolign_restart_subgroups/subgroup_1/pairwise_cigars/ > /private/groups/patenlab/mira/centrolign/guide_tree_testing/method1/chr12_initial_test_nogaps/centrolign_restart_subgroups/subgroup_1/infer_tree/subgroup1_inferred_tree.nwk
```
#### 5. New script to replace subtree

Running on personal computer due to conda issues
```
cd /Users/miramastoras/Desktop/replace_subtrees/

python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/replace_subtree.py \
    -i /Users/miramastoras/Desktop/tree_heatmap_chr12/samples.txt \
    -s subgroup4_inferred_tree.nwk.txt,subgroup0_inferred_tree.nwk.txt,subgroup2_inferred_tree.nwk.txt,subgroup1_inferred_tree.nwk \
    -t KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.nwk.txt \
    -o KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.all_samples.tree_method1.all_subgroups.nwk.txt
```

Re-run pairwise heatmap
```
# all samples
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/replace_subtrees/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.all_samples.tree_method1.all_subgroups.nwk.txt \
        -s /Users/miramastoras/Desktop/tree_heatmap_chr12/samples.txt  \
        -p /Users/miramastoras/Desktop/tree_heatmap_chr12/pairwise_distance.csv \
        -o /Users/miramastoras/Desktop/replace_subtrees/all_subgroups_tree_method1
```
![reinferred](pics/all_subgroups_tree_method1_pairwise_tree_heatmap.png)

#### 6. Re-run centrolign with the newly inferred tree

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
    -S /private/groups/patenlab/mira/centrolign/guide_tree_testing/method1/chr12_initial_test_nogaps/rerun_centrolign_w_new_tree/jobstore/ \
    -T /private/groups/patenlab/mira/centrolign/guide_tree_testing/method1/chr12_initial_test_nogaps/rerun_centrolign_w_new_tree/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.all_samples.tree_method1.all_subgroups.nwk.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/initial_test_no_gaps_chr12.fasta \
    > /private/groups/patenlab/mira/centrolign/guide_tree_testing/method1/chr12_initial_test_nogaps/rerun_centrolign_w_new_tree/initial_test_no_gaps_chr12.tree_method1.centrolign.gfa
```
