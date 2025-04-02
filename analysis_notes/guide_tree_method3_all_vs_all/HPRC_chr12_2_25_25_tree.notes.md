# Testing new guide trees Sasha produced from the HPRC assemblies

First step I needed to re-run the extract_hors wdl to make sure the fasta naming convention matches the HPRC.

https://github.com/miramastoras/centrolign_analysis/tree/main/batch_submissions/extract_hors_HPRC

## Run centrolign directly on the guide trees

Get lists of fasta files containing complete HOR for each chrom
```
cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/initial_test_nogaps

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/
ls | grep hap | while read line ; do
  realpath $line/analysis/extract_hors_HPRC_outputs/*.fasta
done | grep chr12 > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/fasta_list.txt

for f in 100kb 500kb ; do
    cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/initial_test_nogaps_${f}_flanks
    mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps_${f}_flanks/chr12
    ls | grep hap | while read line ; do
    realpath $line/analysis/extract_hors_HPRC_outputs/*.fasta
    done | grep chr12 > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps_${f}_flanks/chr12/fasta_list.txt
  done
```
Get list of samples in nwk tree

```
# get list of sample names as they'd be listed in tree
cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/fasta_list.txt | while read line ; do basename $line | cut -f3 -d"_" ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/fasta_list.all_sample_ids.txt

SAMPLES=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/fasta_list.all_sample_ids.txt
NWK=/private/groups/patenlab/mira/centrolign/annotations/guide_trees/HPRC_chr12_P_Q_mira.2.25.25.rnj.nwk
    while IFS= read -r pattern; do
      if grep -q "$pattern" $NWK; then
        echo "$pattern"
        fi
    done < $SAMPLES > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/sample_ids.in_nwk.txt

# combine fastas
cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/sample_ids.in_nwk.txt | while read line ; do
    grep $line /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/fasta_list.txt
    done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/fasta_list.inside_nwk.txt

# combine fastas
cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/fasta_list.inside_nwk.txt | while read line ; do cat $line ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/initial_test_nogaps_HPRC_labels_chr12.fasta
samtools faidx /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/initial_test_nogaps_HPRC_labels_chr12.fasta
```
Prepare combined fastas for HORs + flanks. Using same list of samples in newick trees.
```
# get list of fasta files
for f in 100kb 500kb ; do
  cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/sample_ids.in_nwk.txt | while read line ; do
    grep $line /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps_${f}_flanks/chr12/fasta_list.txt
    done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps_${f}_flanks/chr12/fasta_list.inside_nwk.txt
done

# combine fastas
for f in 100kb 500kb ; do
    cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps_${f}_flanks/chr12/fasta_list.inside_nwk.txt | while read line ; do cat $line ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps_${f}_flanks/chr12/initial_test_nogaps_${f}_flanks_HPRC_labels_chr12.fasta
    samtools faidx /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps_${f}_flanks/chr12/initial_test_nogaps_${f}_flanks_HPRC_labels_chr12.fasta
  done
```
Run centrolign with new upgma and new nj trees

Convert them to format 5 in ete3 which is accepted by centrolign
```py
from ete3 import Tree

tree = Tree("/Users/miramastoras/Desktop/HPRC_chr12_P_Q_mira.2.25.25.rnj.nwk", format=1)
tree.write(outfile="/Users/miramastoras/Desktop/HPRC_chr12_P_Q_mira.2.25.25.rnj.format5.nwk", format=5)

tree = Tree("/Users/miramastoras/Desktop/HPRC_chr12_P_Q_mira.2.25.25.upgma.nwk", format=1)
tree.write(outfile="/Users/miramastoras/Desktop/HPRC_chr12_P_Q_mira.2.25.25.upgma.format5.nwk", format=5)
```

```
#!/bin/bash
#SBATCH --job-name=centrolign_chr12_2_25_rnj
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=700gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=centrolign_%x.%j.log
#SBATCH --time=7-00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/upgma_tree/jobstore/ \
    -T /private/groups/patenlab/mira/centrolign/annotations/guide_trees/HPRC_chr12_P_Q_mira.2.25.25.upgma.format5.nwk \
    -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/upgma_tree/pairwise_cigars/pairwise_cigar \
    -R \
    --threads 32 \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/initial_test_nogaps_HPRC_labels_chr12.fasta \
    > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/upgma_tree/initial_test_nogaps_HPRC_labels_chr12_2_25_25.upgma.centrolign.gfa

time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/rnj_tree/jobstore/ \
    -T /private/groups/patenlab/mira/centrolign/annotations/guide_trees/HPRC_chr12_P_Q_mira.2.25.25.rnj.format5.nwk \
    -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/rnj_tree/pairwise_cigars/pairwise_cigar \
    -R \
    --threads 32 \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/initial_test_nogaps_HPRC_labels_chr12.fasta \
    > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/rnj_tree/initial_test_nogaps_HPRC_labels_chr12_2_25_25.rnj.centrolign.gfa
```

RUn centrolign with 100 kb flanks
```
#!/bin/bash
#SBATCH --job-name=centrolign_chr12_upgma_100kb
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=700gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign_%x.%j.log
#SBATCH --time=7-00:00

array_job_7828600_task_9.log
time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps_100kb_flanks/chr12/upgma_tree/jobstore/ \
    -T /private/groups/patenlab/mira/centrolign/annotations/guide_trees/HPRC_chr12_P_Q_mira.2.25.25.upgma.format5.nwk \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps_100kb_flanks/chr12/initial_test_nogaps_100kb_flanks_HPRC_labels_chr12.fasta \
    > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps_100kb_flanks/chr12/upgma_tree/initial_test_nogaps_100kb_flanks_HPRC_labels_chr12_2_25_25.upgma.centrolign.gfa

```
## Combine distances from new guide trees with centrolign all pairs alignments, rerun centrolign

Rerun centrolign all pairs pairwise with the new HOR fastas
```
#!/bin/bash
# Slurm script to use centrolign as a pairwise aligner on all pairs of sequences
# from an input directory

#SBATCH --job-name=pairwise-centrolign_hprc
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=56gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=[1-7750]%75
#SBATCH --output=logs/array_job_%A_task_%a.log
#SBATCH --time=1:00:00

date
hostname
pwd

#FLANKS=""
FLANKS=_100kb_flanks
FLANKS=_500kb_flanks
CHROMDIR=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps${FLANKS}/chr12/all_pairs_pairwise/
FASTADIR=/private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/initial_test_nogaps${FLANKS}/
WORKDIR=$CHROMDIR/work/
OUTDIR=$CHROMDIR/pairwise_cigar/

mkdir -p $OUTDIR
mkdir -p $WORKDIR
mkdir -p $FASTADIR

cd $WORKDIR

SAMPLE1=$(awk "NR==$SLURM_ARRAY_TASK_ID" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12_all_pairs/all_pairs_samples.txt | cut -f1)
SAMPLE2=$(awk "NR==$SLURM_ARRAY_TASK_ID" /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12_all_pairs/all_pairs_samples.txt | cut -f2)
echo "sample 1:" $SAMPLE1
echo "sample 2:" $SAMPLE2
echo "out:" $OUTDIR

S1_hap=`echo $SAMPLE1 | cut -f2 -d"."`
S2_hap=`echo $SAMPLE2 | cut -f2 -d"."`

S1_ID=`echo $SAMPLE1 | cut -f1 -d"."`
S2_ID=`echo $SAMPLE2 | cut -f1 -d"."`

if [ ${S1_hap} == 1 ];
then
    S1_HPRC_hap="hap1"
else
    S1_HPRC_hap="hap2"
fi

if [ ${S2_hap} == 1 ];
then
    S2_HPRC_hap="hap1"
else
    S2_HPRC_hap="hap2"
fi

FASTA1=${FASTADIR}${S1_ID}_${S1_HPRC_hap}/analysis/extract_hors_HPRC_outputs/${S1_ID}_${S1_HPRC_hap}_${SAMPLE1}_chr12_hor_array.fasta

FASTA2=${FASTADIR}${S2_ID}_${S2_HPRC_hap}/analysis/extract_hors_HPRC_outputs/${S2_ID}_${S2_HPRC_hap}_${SAMPLE2}_chr12_hor_array.fasta

echo "fasta 1:" $FASTA1
echo "fasta 2:" $FASTA2
TEMP_FASTA=${WORKDIR}/${SAMPLE1}_${SAMPLE2}.fa
cat $FASTA1 $FASTA2 > $TEMP_FASTA
time /private/home/mmastora/progs/centrolign/build/centrolign -v 3 --skip-calibration $TEMP_FASTA > $OUTDIR/pairwise_cigar_${SAMPLE1}_${SAMPLE2}.txt
rm $TEMP_FASTA
```
Got failures for a few samples.

## Testing just trio samples using previous all vs all pairwise runs

infer tree using combination of new SNP distances from Sasha and centrolign all pairs runs for trio samples

swap naming convention
```
sed 's/\.1:/MIRA/g' /Users/miramastoras/Desktop/centrolign_all_pairs_flanks/initial_test_nogaps_100kb_flanks_pairwise_distance.csv | sed 's/\.2:/.1:/g' | sed 's/MIRA/.2:/g' > /Users/miramastoras/Desktop/centrolign_all_pairs_flanks/initial_test_nogaps_100kb_flanks_pairwise_distance.HPRC_convention.csv

# convert to tab separated, removed header line
sed 's/ /\t/g' /Users/miramastoras/Desktop/Sasha_HPRC_trees_2_25_25/HPRC_chr12_P_Q_mira.2.25.25.m > /Users/miramastoras/Desktop/Sasha_HPRC_trees_2_25_25/HPRC_chr12_P_Q_mira.2.25.25.m.tsv
```
```
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/combine_HOR_flank_dist.py \
    -c /Users/miramastoras/Desktop/centrolign_all_pairs_flanks/initial_test_nogaps_100kb_flanks_pairwise_distance.HPRC_convention.csv \
    -f /Users/miramastoras/Desktop/Sasha_HPRC_trees_2_25_25/HPRC_chr12_P_Q_mira.2.25.25.m.tsv \
    -o /Users/miramastoras/Desktop/all_pairs_flanks_tree/100kb_flanks_all_pairs_with_Sasha_dist_HPRC_2_25_25_weighted_sum
```
getting errors with the script because the matrix is differently formatted

### Plot all pairs from centrolign against new HPRC tree
```
# swap naming convention to match Sasha's tree
sed 's/\.1/MIRA/g' /Users/miramastoras/Desktop/centrolign_all_pairs/pairwise_distance.csv | sed 's/\.2/.1/g' | sed 's/MIRA/.2/g' > /Users/miramastoras/Desktop/centrolign_all_pairs/pairwise_distance.HPRC_convention.csv

sed 's/\.1:/MIRA/g' /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.txt | sed 's/\.2:/.1:/g' | sed 's/MIRA/.2:/g' > /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.HPRC_convention.txt
```
Plot tree heatmap
```
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
      -t /Users/miramastoras/Desktop/Sasha_HPRC_trees_2_25_25/HPRC_chr12_P_Q_mira.2.25.25.upgma.format5.nwk \
      -s /Users/miramastoras/Desktop/Sasha_HPRC_trees_2_25_25/samples.txt \
      -p /Users/miramastoras/Desktop/centrolign_all_pairs/pairwise_distance.HPRC_convention.csv \
      -o /Users/miramastoras/Desktop/all_pairs_vs_cenhap_trees/HPRC_chr12_P_Q_mira.2.25.25.upgma_vs_all_pairs_centrolign_
```
Getting this error again
```
File "/Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py", line 95, in tree_to_linkage_matrix
  id2, n_leaves2 = child2.comment
  ^^^^^^^^^^^^^^
TypeError: cannot unpack non-iterable NoneType object
```
Scale branch lengths between 0 and 1
```py
from ete3 import Tree

# Load your tree (replace the string with your actual Newick format)
#tree = Tree("/Users/miramastoras/Desktop/Sasha_HPRC_trees_2_25_25/HPRC_chr12_P_Q_mira.2.25.25.upgma.nwk")
tree = Tree("/Users/miramastoras/Desktop/Sasha_HPRC_trees_2_25_25/HPRC_chr12_P_Q_mira.2.25.25.rnj.format5.nwk")

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
tree.write(outfile="/Users/miramastoras/Desktop/Sasha_HPRC_trees_2_25_25/HPRC_chr12_P_Q_mira.2.25.25.rnj.rescaled_0_1.nwk", format=1)
```
Plot tree heatmap
```
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
      -t /Users/miramastoras/Desktop/Sasha_HPRC_trees_2_25_25/HPRC_chr12_P_Q_mira.2.25.25.upgma.rescaled_0_1.nwk \
      -s /Users/miramastoras/Desktop/Sasha_HPRC_trees_2_25_25/chr12_samples_HPRC_convention.txt \
      -p /Users/miramastoras/Desktop/centrolign_all_pairs/pairwise_distance.HPRC_convention.csv \
      -o /Users/miramastoras/Desktop/all_pairs_vs_cenhap_trees/HPRC_chr12_P_Q_mira.2.25.25.upgma_vs_all_pairs_centrolign_

#
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
      -t /Users/miramastoras/Desktop/Sasha_HPRC_trees_2_25_25/HPRC_chr12_P_Q_mira.2.25.25.rnj.rescaled_0_1.nwk \
      -s /Users/miramastoras/Desktop/Sasha_HPRC_trees_2_25_25/chr12_samples_HPRC_convention.txt \
      -p /Users/miramastoras/Desktop/centrolign_all_pairs/pairwise_distance.HPRC_convention.csv \
      -o /Users/miramastoras/Desktop/all_pairs_vs_cenhap_trees/HPRC_chr12_P_Q_mira.2.25.25.rnj_vs_all_pairs_centrolign_
```

```
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
      -t /Users/miramastoras/Desktop/tree_heatmap/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.nwk.txt \
      -s /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.txt \
      -p /Users/miramastoras/Desktop/centrolign_all_pairs/pairwise_distance.csv \
      -o /Users/miramastoras/Desktop/all_pairs_vs_cenhap_trees/KGP4_TRIOS.upgma_vs_all_pairs_centrolign_
```

Compare sasha's old upgma tree and new one
Swap to old sample ids
```
sed 's/\.1:/MIRA/g' /Users/miramastoras/Desktop/replace_subtrees/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.chr12_hprc_samples.nwk.txt | sed 's/\.2:/.1:/g' | sed 's/MIRA/.2:/g' > /Users/miramastoras/Desktop/replace_subtrees/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.chr12_hprc_samples.HPRC_naming.nwk.txt
```

```
install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
library(phytools)
library(ape)

# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/replace_subtrees/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.chr12_hprc_samples.HPRC_naming.nwk.txt")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/Sasha_HPRC_trees_2_25_25/HPRC_chr12_P_Q_mira.2.25.25.upgma.nwk")

samples <- readLines("/Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.HPRC_naming.txt")

obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/Sasha_HPRC_trees_2_25_25/KGP4_TRIOS_MAC5_chr12_vs_HPRC_chr12_P_Q_mira.2.25.25.upgma.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()
```
Plot centrolign MSA distances against the new RNJ and UPGMA trees

convert cigars to distance
```
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/cigar_to_distance.py /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/upgma_tree/pairwise_cigars/

python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/cigar_to_distance.py /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_chr12_2_25_25_tree_initial_test_nogaps/chr12/rnj_tree/pairwise_cigars/
```

Plot tree_heatmap
```
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
      -t /Users/miramastoras/Desktop/Sasha_HPRC_trees_2_25_25/HPRC_chr12_P_Q_mira.2.25.25.upgma.rescaled_0_1.nwk \
      -s /Users/miramastoras/Desktop/centrolign_chr12_new_trees_2_25_25/sample_ids.in_nwk.new_tree_upgma.txt \
      -p /Users/miramastoras/Desktop/centrolign_chr12_new_trees_2_25_25/upgma_pairwise_distance.csv \
      -o /Users/miramastoras/Desktop/centrolign_chr12_new_trees_2_25_25/centrolign_MSA_induced_pairwise_upgma_tree

s
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
      -t /Users/miramastoras/Desktop/Sasha_HPRC_trees_2_25_25/HPRC_chr12_P_Q_mira.2.25.25.rnj.rescaled_0_1.nwk \
      -s /Users/miramastoras/Desktop/centrolign_chr12_new_trees_2_25_25/sample_ids.in_nwk.new_tree_upgma.txt \
      -p /Users/miramastoras/Desktop/centrolign_chr12_new_trees_2_25_25/rnj_pairwise_distance.csv \
      -o /Users/miramastoras/Desktop/centrolign_chr12_new_trees_2_25_25/centrolign_MSA_induced_pairwise_rnj_tree
```
