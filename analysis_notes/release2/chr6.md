## Centrolign analysis of chr 6 for release 2

### 0. Prepare data

#### Cenhap distances

Cenhap guide trees and distance matrix provided by Chuck and Sasha Langley on 4/11/2025. These 7 samples are removed due to large deletions:
```
NA20827.1
HG00642.2
HG00344.2
HG00639.1
HG02622.1
HG02145.1
NA20799.2
```
distance matrix from p and q
```sh
/private/groups/patenlab/mira/centrolign/annotations/guide_trees/HPRC_chr6_58200000_58286706_61058390_61121735_het68_m_mira_dgp_rnj.m
```
Reformat matrix file to enable reading into python script
```sh
# convert spaces to comma, remove trailing commas
sed 's/ /,/g' HPRC_chr6_58200000_58286706_61058390_61121735_het68_m_mira_dgp_rnj.m | sed 's/,$//' | sed 's/\t/,/' | tail -n +2 > matrix_no_header.csv

# convert row names to column names
N=184
cut -f1 -d"," matrix_no_header.csv | sed 's/$/,/g' | tr -d '\n' | cut -f1-$N -d"," | sed 's/^/0,/g' > rowname_to_colname.csv

cat rowname_to_colname.csv matrix_no_header.csv > HPRC_chr6_58200000_58286706_61058390_61121735_het68_m_mira_dgp_rnj.reformatted.m

rm matrix_no_header.csv rowname_to_colname.csv
```
#### Centrolign all pairs distances

Convert cigar to distance matrix
```sh
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/cigar_to_distance.py /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr6/pairwise_cigar/

# header required for next script
echo "sample1,sample2,distance" > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr6/pairwise_distance.csv

# remove sample not in tree
grep -v "NA20827.1" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr6/pairwise_cigar/pairwise_distance.csv | grep -v "HG00642.2" | grep -v "HG00344.2" | grep -v "HG00639.1" | grep -v "HG02622.1" | grep -v "HG02145.1" | grep -v "NA20799.2" >> /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr6/pairwise_distance.csv
```
### 1. Create refined tree combining centrolign all pairs distances and cenhap distances

Run script to combine distances with weighted sum and compute NJ tree
```sh
conda activate skbio

python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/combine_HOR_flank_dist.py \
    -c /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr6/pairwise_distance.csv \
    -f /private/groups/patenlab/mira/centrolign/annotations/guide_trees/HPRC_chr6_58200000_58286706_61058390_61121735_het68_m_mira_dgp_rnj.reformatted.m \
    -o /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_weighted_sum/HPRC_r2_chr6_cenhap_20250414_centrolign_all_pairs
```
Convert tree to ete3 format 5 which is accepted by centrolign

```py
from ete3 import Tree

tree = Tree("/private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_weighted_sum/HPRC_r2_chr6_cenhap_20250414_centrolign_all_pairs_HOR_flank_dist_weighted.nwk", format=1)
tree.write(outfile="/private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_weighted_sum/HPRC_r2_chr6_cenhap_20250414_centrolign_all_pairs_HOR_flank_dist_weighted.format5.nwk", format=5)
```
Remove 7 samples from fasta file
```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2

cat /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/HPRC_release2_contiguous_HORs_chr6.fasta_list.txt | grep -v "HG00642.2" | grep -v "HG00344.2" | grep -v "HG00639.1" | grep -v "HG02622.1" | grep -v "HG02145.1" | grep -v "NA20799.2" | grep -v "NA20827.1" | while read line ; do cat $line ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/centrolign_fastas/HPRC_release2_HORs_chr6.excl_7_del_samples.fasta
samtools faidx /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/centrolign_fastas/HPRC_release2_HORs_chr6.excl_7_del_samples.fasta
```
Run centrolign on refined tree
```sh
#!/bin/bash
#SBATCH --job-name=centrolign_release2_chr6_refined_tree
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
    -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/MSA_refined_tree/chr6/jobstore/ \
    -T /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_weighted_sum/HPRC_r2_chr6_cenhap_20250414_centrolign_all_pairs_HOR_flank_dist_weighted.format5.nwk \
    /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/centrolign_fastas/HPRC_release2_HORs_chr6.excl_7_del_samples.fasta \
    > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/MSA_refined_tree/chr6/HPRC_r2.chr6.refined_tree.centrolign.gfa
```

Infer a tree from centrolign all pairs alignments of HOR
```sh
conda activate skbio

python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/data_exploration/infer_tree.py /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr6/pairwise_cigar/ > /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_all_pairs/chr6_r2_centrolign_all_pairs_nj_tree.nwk
```

Reformat for centrolign
```py
from ete3 import Tree

tree = Tree("/private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_all_pairs/chr6_r2_centrolign_all_pairs_nj_tree.nwk", format=1)
tree.write(outfile="/private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_all_pairs/chr6_r2_centrolign_all_pairs_nj_tree.format5.nwk", format=5)
```
Run centrolign
```sh
#!/bin/bash
#SBATCH --job-name=centrolign_release2_chr6_allpairs_NJ_HOR_tree
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
    -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/MSA_HOR_all_pairs_NJ_tree/chr6/jobstore/ \
    -T /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_all_pairs/chr6_r2_centrolign_all_pairs_nj_tree.format5.nwk \
    /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/centrolign_fastas/HPRC_release2_HORs_chr6.excl_7_del_samples.fasta \
    > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/MSA_HOR_all_pairs_NJ_tree/chr6/HPRC_r2.chr6.allpairs_HOR_NJ.centrolign.gfa
```

Plot tanglegram between refined tree and HOR tree:
```
install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
library(phytools)
library(ape)

# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/chr6_r2_centrolign_all_pairs_nj_tree.format5.nwk")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/HPRC_r2_chr6_cenhap_20250414_centrolign_all_pairs_HOR_flank_dist_weighted.format5.nwk")

obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/chr6_r2_all_pairs_vs_refined_tree.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()
```
