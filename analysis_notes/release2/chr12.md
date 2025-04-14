## Centrolign analysis of chr 12 for release 2

- Using chr12 as a first pass to implement our analysis plan, before automating it for the other chroms.

### 0. Prepare data

#### Cenhap distances

Cenhap guide trees and distance matrix provided by Chuck and Sasha Langley on 4/2/2025, containing all the samples with contiguous HORs from release 2, except HG00741.1 because it has a large deletion(?) relative to chm13 on the q side of the array

```sh
/private/groups/patenlab/mira/centrolign/annotations/guide_trees/HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj.m
/private/groups/patenlab/mira/centrolign/annotations/guide_trees/HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj.nwk
/private/groups/patenlab/mira/centrolign/annotations/guide_trees/HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj_upgma.nwk
```

Reformat matrix file to enable reading into python script
```
# convert spaces to comma, remove trailing commas
sed 's/ /,/g' HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj.m | sed 's/,$//' | sed 's/\t//' | tail -n +2 > matrix_no_header.csv

# convert row names to column names
N=384
cut -f1 -d"," matrix_no_header.csv | sed 's/$/,/g' | tr -d '\n' | cut -f1-$N -d"," | sed 's/^/0,/g' > rowname_to_colname.csv

cat rowname_to_colname.csv matrix_no_header.csv > HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj.reformatted.m

rm matrix_no_header.csv rowname_to_colname.csv
```
#### Centrolign all pairs distances

Convert cigar to distance matrix
```
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/cigar_to_distance.py /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr12/pairwise_cigar/

# header required for next script
echo "sample1,sample2,distance" > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr12/pairwise_distance_excl_HG00741.1.csv

# remove sample not in tree
grep -v "HG00741.1" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr12/pairwise_distance.csv >> /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr12/pairwise_distance_excl_HG00741.1.csv
```

### 1. Create refined tree combining centrolign all pairs distances and cenhap distances

Run script to combine distances with weighted sum and compute NJ tree
```
conda activate skbio

python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/combine_HOR_flank_dist.py \
    -c /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr12/pairwise_distance_excl_HG00741.1.csv \
    -f /private/groups/patenlab/mira/centrolign/annotations/guide_trees/HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj.reformatted.m \
    -o /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_weighted_sum/HPRC_r2_chr12_cenhap_20250402_centrolign_all_pairs
```
Convert tree to ete3 format 5 which is accepted by centrolign

```py
from ete3 import Tree

tree = Tree("/private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_weighted_sum/HPRC_r2_chr12_cenhap_20250402_centrolign_all_pairs_HOR_flank_dist_weighted.nwk", format=1)
tree.write(outfile="/private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_weighted_sum/HPRC_r2_chr12_cenhap_20250402_centrolign_all_pairs_HOR_flank_dist_weighted.format5.nwk", format=5)
```

Remove HG00741.1 from fasta file
```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2

cat /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/HPRC_release2_contiguous_HORs_chr12.fasta_list.txt | grep -v "HG00741.1" | while read line ; do cat $line ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/centrolign_fastas/HPRC_release2_HORs_chr12.excl_HG00741.1.fasta
samtools faidx /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/centrolign_fastas/HPRC_release2_HORs_chr12.excl_HG00741.1.fasta
```
Run centrolign on refined tree
```sh
#!/bin/bash
#SBATCH --job-name=centrolign_release2_chr12_refined_tree
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
    -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/MSA_refined_tree/chr12/jobstore/ \
    -T /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_weighted_sum/HPRC_r2_chr12_cenhap_20250402_centrolign_all_pairs_HOR_flank_dist_weighted.format5.nwk \
    /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/centrolign_fastas/HPRC_release2_HORs_chr12.excl_HG00741.1.fasta \
    > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/MSA_refined_tree/chr12/HPRC_r2.chr12.refined_tree.centrolign.gfa
```
Restart with pairwise cigars multithreaded:
```
-A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/upgma_tree/pairwise_cigars/pairwise_cigar \
-R \
--threads 32 \
```

Run centrolign MSA on NJ cenhap tree only

Fix weird spacing issue
```sh
sed "s/ '//g" HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj.nwk | sed "s/'//g"
```
Reformat for centrolign
```py
from ete3 import Tree

tree = Tree("/Users/miramastoras/Desktop/HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj.fix.nwk", format=1)
tree.write(outfile="/Users/miramastoras/Desktop/HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj.format5.nwk", format=5)
```

Run centrolign
```sh
#!/bin/bash
#SBATCH --job-name=centrolign_release2_chr12_NJ_cenhap_tree
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
    -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/MSA_NJ_cenhap_tree/chr12/jobstore/ \
    -T /private/groups/patenlab/mira/centrolign/annotations/guide_trees/HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj.format5.nwk \
    /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/centrolign_fastas/HPRC_release2_HORs_chr12.excl_HG00741.1.fasta \
    > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/MSA_NJ_cenhap_tree/chr12/HPRC_r2.chr12.cenhap_NJ.centrolign.gfa
```

Infer a tree from centrolign all pairs alignments of HOR
```sh
conda activate skbio

python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/data_exploration/infer_tree.py /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr12/pairwise_cigar/ > /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_all_pairs/chr12_r2_centrolign_all_pairs_nj_tree.nwk
```

Reformat for centrolign
```py
from ete3 import Tree

tree = Tree("/Users/miramastoras/Desktop/chr12_r2_centrolign_all_pairs_nj_tree.nwk", format=1)
tree.write(outfile="/Users/miramastoras/Desktop/chr12_r2_centrolign_all_pairs_nj_tree.format5.nwk", format=5)
```
Run centrolign
```sh
#!/bin/bash
#SBATCH --job-name=centrolign_release2_chr12_allpairs_NJ_HOR_tree
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
    -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/MSA_HOR_all_pairs_NJ_tree/chr12/jobstore/ \
    -T /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_all_pairs/chr12_r2_centrolign_all_pairs_nj_tree.format5.nwk \
    /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/centrolign_fastas/HPRC_release2_HORs_chr12.excl_HG00741.1.fasta \
    > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/MSA_HOR_all_pairs_NJ_tree/chr12/HPRC_r2.chr12.allpairs_HOR_NJ.centrolign.gfa
```

Run centrolign on UPGMA cenhap tree as well

```py
# Reformat for centrolign
from ete3 import Tree

tree = Tree("/Users/miramastoras/Desktop/HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj_upgma.nwk", format=1)
tree.write(outfile="/Users/miramastoras/Desktop/HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj_upgma.nwk.format5.nwk", format=5)
```
Submit centrolign
```sh
#!/bin/bash
#SBATCH --job-name=centrolign_release2_chr12_upgma_cenhap_tree
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
    -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/MSA_UPGMA_cenhap_tree/chr12/jobstore/ \
    -T /private/groups/patenlab/mira/centrolign/annotations/guide_trees/HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj_upgma.nwk.format5.nwk \
    /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/centrolign_fastas/HPRC_release2_HORs_chr12.excl_HG00741.1.fasta \
    > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/MSA_UPGMA_cenhap_tree/chr12/HPRC_r2.chr12.upgma_cenhap_tree.centrolign.gfa
```
## 2. Tree evaluation

### 2.1 Plot pairwise heatmap
NJ cenhap tree

Scale branch lengths to be between 0 and 1
```py
from ete3 import Tree

tree = Tree("/Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj.format5.nwk")

# A small epsilon to avoid zero branch lengths
epsilon = 1e-6

# Get all the branch lengths
branch_lengths = [n.dist for n in tree.traverse()]

# Shift the branch lengths by adding epsilon to avoid 0, and then rescale
min_length = min(branch_lengths)
max_length = max(branch_lengths)

# If all branch lengths are the same, scaling won't work. In that case, just return 1 for all.
if min_length == max_length:
    for node in tree.traverse():
        node.dist = 1
else:
    # Scale the branch lengths to be between epsilon and 1
    def scale_branch_length(length):
        # Add epsilon to avoid zero and rescale to the range (0, 1)
        return (length - min_length + epsilon) / (max_length - min_length + epsilon)

    # Apply the scaling function to each branch length
    for node in tree.traverse():
        node.dist = scale_branch_length(node.dist)

# Check the tree with scaled branch lengths
# tree.write(outfile="/Users/miramastoras/Desktop/Sasha_HPRC_trees_2_25_25/HPRC_chr12_P_Q_mira.2.25.25.upgma.rescaled_0_1.nwk", format=1)
tree.write(outfile="/Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj.format5.rescaled_0_1.nwk", format=5)

```

```sh
conda activate tree_python

python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj.format5.rescaled_0_1.nwk \
        -s /Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_release2_contiguous_HORs_chr12.excl_HG00741.1.top130.txt \
        -p /Users/miramastoras/Desktop/chr12_r2_tree_comparison/pairwise_distance_excl_HG00741.1.csv \
        -m "Centrolign all pairs" \
        -n "NJ Cenhap Tree" \
        -d "NJ Distances" \
        -o /Users/miramastoras/Desktop/chr12_r2_tree_comparison/chr12_centrolign_all_pairs.cenhap_nj_tree_

# error TypeError: cannot unpack non-iterable NoneType object solved by scaling tree branch lengths between 0 and 1

```
UPGMA cenhap tree

Scale branch lengths to be between 0 and 1
```py
from ete3 import Tree

tree = Tree("/Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj_upgma.nwk")

# A small epsilon to avoid zero branch lengths
epsilon = 1e-6

# Get all the branch lengths
branch_lengths = [n.dist for n in tree.traverse()]

# Shift the branch lengths by adding epsilon to avoid 0, and then rescale
min_length = min(branch_lengths)
max_length = max(branch_lengths)

# If all branch lengths are the same, scaling won't work. In that case, just return 1 for all.
if min_length == max_length:
    for node in tree.traverse():
        node.dist = 1
else:
    # Scale the branch lengths to be between epsilon and 1
    def scale_branch_length(length):
        # Add epsilon to avoid zero and rescale to the range (0, 1)
        return (length - min_length + epsilon) / (max_length - min_length + epsilon)

    # Apply the scaling function to each branch length
    for node in tree.traverse():
        node.dist = scale_branch_length(node.dist)

# Check the tree with scaled branch lengths
# tree.write(outfile="/Users/miramastoras/Desktop/Sasha_HPRC_trees_2_25_25/HPRC_chr12_P_Q_mira.2.25.25.upgma.rescaled_0_1.nwk", format=1)
tree.write(outfile="/Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj_upgma.rescaled_0_1.nwk", format=5)

```


```sh
conda activate tree_python

python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj_upgma.rescaled_0_1.nwk \
        -s /Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_release2_contiguous_HORs_chr12.excl_HG00741.1.top130.txt \
        -p /Users/miramastoras/Desktop/chr12_r2_tree_comparison/pairwise_distance_excl_HG00741.1.csv \
        -m "Centrolign all pairs" \
        -n "UPGMA Cenhap Tree" \
        -d "UPGMA Distances" \
        -o /Users/miramastoras/Desktop/chr12_r2_tree_comparison/chr12_centrolign_all_pairs.cenhap_upgma_tree_

# error TypeError: cannot unpack non-iterable NoneType object

```
HOR tree
```sh
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/chr12_r2_tree_comparison/chr12_r2_centrolign_all_pairs_nj_tree.format5.nwk \
        -s /Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_release2_contiguous_HORs_chr12.excl_HG00741.1.top130.txt \
        -p /Users/miramastoras/Desktop/chr12_r2_tree_comparison/pairwise_distance_excl_HG00741.1.csv \
        -m "Centrolign all pairs distances" \
        -n "Centrolign all pairs HOR NJ Tree" \
        -d "All pairs Distances" \
        -o /Users/miramastoras/Desktop/chr12_r2_tree_comparison/chr12_centrolign_all_pairs.HOR_NJ_tree_

# ValueError: All values in the dash list must be non-negative

```

Refined Tree
```sh
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_r2_chr12_cenhap_20250402_centrolign_all_pairs_HOR_flank_dist_weighted.format5.nwk \
        -s /Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_release2_contiguous_HORs_chr12.excl_HG00741.1.top200.txt \
        -p /Users/miramastoras/Desktop/chr12_r2_tree_comparison/pairwise_distance_excl_HG00741.1.csv \
        -m "Centrolign all pairs distances" \
        -n "Refined tree (HOR + cenhap flanks)" \
        -d "Tree Distances" \
        -o /Users/miramastoras/Desktop/chr12_r2_tree_comparison/chr12_centrolign_refined_tree
```
### 2.2 Run linear regression on all the trees

Cenhap NJ tree

```sh
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/regression_tree_eval.py \
    -t /Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj.format5.rescaled_0_1.nwk \
    -s /Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_release2_contiguous_HORs_chr12.excl_HG00741.1.txt \
    -d /Users/miramastoras/Desktop/chr12_r2_tree_comparison/pairwise_distance_excl_HG00741.1.csv \
    -f col \
    -o /Users/miramastoras/Desktop/chr12_r2_tree_comparison/cenhap_NJ_tree_centrolign_all_pairs_

# Results:

Number of internal nodes: 382
Number of leaves: 384
Root Mean Squared Error: 0.10367261687517498
```

Cenhap UPGMA tree
```sh
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/regression_tree_eval.py \
    -t /Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj_upgma.rescaled_0_1.nwk \
    -s /Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_release2_contiguous_HORs_chr12.excl_HG00741.1.txt \
    -d /Users/miramastoras/Desktop/chr12_r2_tree_comparison/pairwise_distance_excl_HG00741.1.csv \
    -f col \
    -o /Users/miramastoras/Desktop/chr12_r2_tree_comparison/cenhap_upgma_tree_centrolign_all_pairs_

# Results
Number of internal nodes: 383
Number of leaves: 384
Root Mean Squared Error: 0.033244717502594806
```
HOR only tree
```sh
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/regression_tree_eval.py \
    -t /Users/miramastoras/Desktop/chr12_r2_tree_comparison/chr12_r2_centrolign_all_pairs_nj_tree.format5.nwk \
    -s /Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_release2_contiguous_HORs_chr12.excl_HG00741.1.txt \
    -d /Users/miramastoras/Desktop/chr12_r2_tree_comparison/pairwise_distance_excl_HG00741.1.csv \
    -f col \
    -o /Users/miramastoras/Desktop/chr12_r2_tree_comparison/HOR_NJ_tree_centrolign_all_pairs_

# Results
Number of internal nodes: 383
Number of leaves: 384
Root Mean Squared Error: 0.02714508479152619
```
Refined Tree
```sh
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/regression_tree_eval.py \
    -t /Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_r2_chr12_cenhap_20250402_centrolign_all_pairs_HOR_flank_dist_weighted.format5.nwk \
    -s /Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_release2_contiguous_HORs_chr12.excl_HG00741.1.txt \
    -d /Users/miramastoras/Desktop/chr12_r2_tree_comparison/pairwise_distance_excl_HG00741.1.csv \
    -f col \
    -o /Users/miramastoras/Desktop/chr12_r2_tree_comparison/refined_tree_centrolign_all_pairs_

# Results
Number of internal nodes: 383
Number of leaves: 384
Root Mean Squared Error: 0.0275442624524752
```
#### Plot the residuals

```sh
conda activate tree_python

python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj.format5.rescaled_0_1.nwk \
        -s /Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_release2_contiguous_HORs_chr12.excl_HG00741.1.top130.txt \
        -p /Users/miramastoras/Desktop/chr12_r2_tree_comparison/cenhap_NJ_tree_centrolign_all_pairs_residuals_pairwise.csv \
        -m "Centrolign all pairs residuals" \
        -n "NJ Cenhap Tree" \
        -d "NJ Distances" \
        -o /Users/miramastoras/Desktop/chr12_r2_tree_comparison/chr12_centrolign_all_pairs_residuals.cenhap_nj_tree_

# error TypeError: cannot unpack non-iterable NoneType object solved by scaling tree branch lengths between 0 and 1

```
UPGMA cenhap tree

```sh
conda activate tree_python

python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj_upgma.rescaled_0_1.nwk \
        -s /Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_release2_contiguous_HORs_chr12.excl_HG00741.1.top130.txt \
        -p /Users/miramastoras/Desktop/chr12_r2_tree_comparison/cenhap_upgma_tree_centrolign_all_pairs_residuals_pairwise.csv \
        -m "Centrolign all pairs residuals" \
        -n "UPGMA Cenhap Tree" \
        -d "UPGMA Distances" \
        -o /Users/miramastoras/Desktop/chr12_r2_tree_comparison/chr12_centrolign_all_pairs_residuals.cenhap_upgma_tree_

# error TypeError: cannot unpack non-iterable NoneType object

```
HOR tree
```sh
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/chr12_r2_tree_comparison/chr12_r2_centrolign_all_pairs_nj_tree.format5.nwk \
        -s /Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_release2_contiguous_HORs_chr12.excl_HG00741.1.top130.txt \
        -p /Users/miramastoras/Desktop/chr12_r2_tree_comparison/HOR_NJ_tree_centrolign_all_pairs_residuals_pairwise.csv \
        -m "Centrolign all pairs residuals" \
        -n "Centrolign all pairs HOR NJ Tree" \
        -d "All pairs Distances" \
        -o /Users/miramastoras/Desktop/chr12_r2_tree_comparison/chr12_centrolign_all_pairs_residuals.HOR_NJ_tree_

# ValueError: All values in the dash list must be non-negative

```

Refined Tree
```sh
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_r2_chr12_cenhap_20250402_centrolign_all_pairs_HOR_flank_dist_weighted.format5.nwk \
        -s /Users/miramastoras/Desktop/chr12_r2_tree_comparison/HPRC_release2_contiguous_HORs_chr12.excl_HG00741.1.top130.txt \
        -p /Users/miramastoras/Desktop/chr12_r2_tree_comparison/refined_tree_centrolign_all_pairs_residuals_pairwise.csv \
        -m "Centrolign all pairs residuals" \
        -n "Refined tree (HOR + cenhap flanks)" \
        -d "Tree Distances" \
        -o /Users/miramastoras/Desktop/chr12_r2_tree_comparison/chr12_centrolign_refined_tree_residuals
```
