### Using all vs all pairwise alignments to infer the tree

```
while read -r s1; do
    while read -r s2; do
        [[ "$s1" < "$s2" ]] && echo -e "$s1\t$s2"
    done < /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/samples.txt
done < /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/samples.txt > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12_all_pairs/all_pairs_samples.txt

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12_all_pairs/

mkdir logs
```
slurm script
```
#!/bin/bash
# Slurm script to use centrolign as a pairwise aligner on all pairs of sequences
# from an input directory

#SBATCH --job-name=pairwise-centrolign
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=56gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=[1-7750]%150
#SBATCH --output=logs/array_job_%A_task_%a.log
#SBATCH --time=1:00:00

date
hostname
pwd

CHROMDIR=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12_all_pairs
FASTADIR=/private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors/initial_test_nogaps/
WORKDIR=$CHROMDIR/work/
OUTDIR=$CHROMDIR/pairwise_cigar/

mkdir -p $OUTDIR
mkdir -p $WORKDIR
mkdir -p $FASTADIR

cd $WORKDIR

SAMPLE1=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$CHROMDIR"/all_pairs_samples.txt | cut -f1)
SAMPLE2=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$CHROMDIR"/all_pairs_samples.txt | cut -f2)
echo "sample 1:" $SAMPLE1
echo "sample 2:" $SAMPLE2
echo "out:" $OUTDIR

S1_hap=`echo $SAMPLE1 | cut -f2 -d"."`
S2_hap=`echo $SAMPLE2 | cut -f2 -d"."`

S1_ID=`echo $SAMPLE1 | cut -f1 -d"."`
S2_ID=`echo $SAMPLE2 | cut -f1 -d"."`

if [ ${S1_hap} == 1 ];
then
    S1_HPRC_hap="hap2"
else
    S1_HPRC_hap="hap1"
fi

if [ ${S2_hap} == 1 ];
then
    S2_HPRC_hap="hap2"
else
    S2_HPRC_hap="hap1"
fi

FASTA1=${FASTADIR}${S1_ID}_${S1_HPRC_hap}/analysis/extract_hors_outputs/${S1_ID}_${S1_HPRC_hap}_${SAMPLE1}_chr12_hor_array.fasta

FASTA2=${FASTADIR}${S2_ID}_${S2_HPRC_hap}/analysis/extract_hors_outputs/${S2_ID}_${S2_HPRC_hap}_${SAMPLE2}_chr12_hor_array.fasta

echo "fasta 1:" $FASTA1
echo "fasta 2:" $FASTA2
TEMP_FASTA=${WORKDIR}/${SAMPLE1}_${SAMPLE2}.fa
cat $FASTA1 $FASTA2 > $TEMP_FASTA
time /private/home/mmastora/progs/centrolign/build/centrolign -v 3 --skip-calibration $TEMP_FASTA > $OUTDIR/pairwise_cigar_${SAMPLE1}_${SAMPLE2}.txt
rm $TEMP_FASTA
```
Run neighbor joining tree using these pairwise distances
```
conda activate skbio

python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/data_exploration/infer_tree.py /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12_all_pairs/pairwise_cigar/ > /private/groups/patenlab/mira/centrolign/guide_tree_testing/method3_all_vs_all/chr12_centrolign_all_pairs_nj_tree.nwk

# change naming
ls | while read line ; do
new_name=`echo $line | sed 's\cigar.\pairwise_cigar_\g'`
mv ${line} ${new_name}
done
```
convert to correct format
```
from ete3 import Tree

tree = Tree("/Users/miramastoras/Desktop/centrolign_all_pairs/chr12_centrolign_all_pairs_nj_tree.nwk", format=1)
tree.write(outfile="/Users/miramastoras/Desktop/centrolign_all_pairs/chr12_centrolign_all_pairs_nj_tree.format5.nwk.txt", format=5)
```
Plot against pairwise distances
```
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/cigar_to_distance.py /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12_all_pairs/pairwise_cigar/

# on personal computer
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/centrolign_all_pairs/chr12_centrolign_all_pairs_nj_tree.format5.nwk.txt \
        -s /Users/miramastoras/Desktop/tree_heatmap_chr12/samples.txt \
        -p /Users/miramastoras/Desktop/centrolign_all_pairs/pairwise_distance.csv \
        -o /Users/miramastoras/Desktop/centrolign_all_pairs/chr12_centrolign_all_pairs_nj_tree.all_pairs_distances_

/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12_all_pairs/pairwise_cigar
```
### Using jordan's formula & all vs all pairwise distances to infer initial tree

First subset to only trio samples
```
grep "_mat_" /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors/initial_test_nogaps/extract_hors_initial_test_nogaps.csv >/private/groups/patenlab/mira/centrolign/test/samples_trio.txt
grep "_pat_" /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors/initial_test_nogaps/extract_hors_initial_test_nogaps.csv >> /private/groups/patenlab/mira/centrolign/test/samples_trio.txt

cut -f 1 /private/groups/patenlab/mira/centrolign/test/samples_trio.txt | cut -f 1 -d"_" | while read line ; do grep $line /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/fasta_list.all_sample_ids.in_nwk.txt ; done | sort | uniq > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/fasta_list.all_sample_ids.in_nwk.trio_only.txt
```

```
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/combine_HOR_flank_dist.py \
    -c /Users/miramastoras/Desktop/centrolign_all_pairs/pairwise_distance.csv \
    -f /Users/miramastoras/Desktop/distance_compare/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_hap.hD.txt \
    -o /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/chr12
```
Plot pairwise heatmap
```
conda activate tree_python

# combined distances as a spot check
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/chr12_HOR_flank_dist_weighted.nwk \
        -s /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.txt \
        -p /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/chr12_HOR_flank_dist_weighted.txt \
        -o /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/combine_HOR_flank_dist_chr12_all_pairs_

# just centrolign distances
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/chr12_HOR_flank_dist_weighted.nwk \
        -s /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.txt \
        -p /Users/miramastoras/Desktop/centrolign_all_pairs/pairwise_distance.csv \
        -o /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/combine_HOR_flank_dist_chr12_all_pairs_aln
```

### Test expanding out sequence into the flanks for chr12

Slides: https://docs.google.com/presentation/d/1ucEL2VZXHVe3YCSrSSkUY6CO_osXBAWDeyyitBwrHOU/edit#slide=id.g331a6f2c84d_0_61

Testing 100kb, 50kb, 500kb

Get lists of fasta files containing complete HOR for each chrom
```
# 100 kb
cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors

for f in 100kb 50kb 500kb; do
    mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps_${f}_flanks/chr12/
    ls initial_test_nogaps_${f}_flanks/ | grep hap | \
    while read line ; do
        realpath initial_test_nogaps_${f}_flanks/${line}/analysis/extract_hors_outputs/*.fasta
      done | grep chr12 > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps_${f}_flanks/chr12/fasta_list_complete_HORs.txt
  done
```
Get list of samples that are in the nwk tree
```
# get sample ids in cenhap convention
for f in 100kb 50kb 500kb; do  
    cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps_${f}_flanks/chr12/fasta_list_complete_HORs.txt | while read line ; do basename $line | cut -f3 -d"_" ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps_${f}_flanks/chr12/fasta_list.all_sample_ids.txt
  done

# extract those sample ids in newick tree
for f in 100kb 50kb 500kb; do   
    SAMPLES=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps_${f}_flanks/chr12/fasta_list.all_sample_ids.txt
    NWK=/private/groups/patenlab/mira/centrolign/annotations/cenhap_test/KGP4_TRIOS_MAC5_chr12_CPR_no_filter_EHet30_no_PS_PID_PGT_lifted_over.v1.1_new_10_haps_plus.txt
    while IFS= read -r pattern; do
      if grep -q "$pattern" $NWK; then
        echo "$pattern"
        fi
    done < $SAMPLES > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps_${f}_flanks/chr12/sample_list_complete_HORs.in_nwk.txt
  done

# extract full fasta paths for those sample ids in newick tree
for f in 100kb 50kb 500kb; do    
    cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps_${f}_flanks/chr12/sample_list_complete_HORs.in_nwk.txt | while read line ; do
    grep $line /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps_${f}_flanks/chr12/fasta_list_complete_HORs.txt
    done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps_${f}_flanks/chr12/fasta_list_complete_HORs.in_nwk.txt
done

# combine all sequences into one fasta
for f in 100kb 50kb 500kb; do
    cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps_${f}_flanks/chr12/fasta_list_complete_HORs.in_nwk.txt | while read line ; do cat $line ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps_${f}_flanks/chr12/initial_test_nogaps_${f}_chr12.fasta
    samtools faidx /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps_${f}_flanks/chr12/initial_test_nogaps_${f}_chr12.fasta
done
```
Run centrolign all pairs - pairwise aligner

```
for f in 100kb 50kb 500kb; do
    mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps_${f}_flanks/chr12/all_pairs_pairwise/logs
  done
```
Prepare all pairs sample list
```
for f in 100kb 50kb 500kb; do
  while read -r s1; do
      while read -r s2; do
          [[ "$s1" < "$s2" ]] && echo -e "$s1\t$s2"
          done < /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps_${f}_flanks/chr12/sample_list_complete_HORs.in_nwk.txt
      done < /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps_${f}_flanks/chr12/sample_list_complete_HORs.in_nwk.txt > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps_${f}_flanks/chr12/all_pairs_pairwise/all_pairs_samples.txt

```

```
#!/bin/bash
# Slurm script to use centrolign as a pairwise aligner on all pairs of sequences
# from an input directory

#SBATCH --job-name=pairwise-centrolign_50kb
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

FLANKS=50kb
CHROMDIR=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps_${FLANKS}_flanks/chr12/all_pairs_pairwise/
FASTADIR=/private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors/initial_test_nogaps_${FLANKS}_flanks/
WORKDIR=$CHROMDIR/work/
OUTDIR=$CHROMDIR/pairwise_cigar/

mkdir -p $OUTDIR
mkdir -p $WORKDIR
mkdir -p $FASTADIR

cd $WORKDIR

SAMPLE1=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$CHROMDIR"/all_pairs_samples.txt | cut -f1)
SAMPLE2=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$CHROMDIR"/all_pairs_samples.txt | cut -f2)
echo "sample 1:" $SAMPLE1
echo "sample 2:" $SAMPLE2
echo "out:" $OUTDIR

S1_hap=`echo $SAMPLE1 | cut -f2 -d"."`
S2_hap=`echo $SAMPLE2 | cut -f2 -d"."`

S1_ID=`echo $SAMPLE1 | cut -f1 -d"."`
S2_ID=`echo $SAMPLE2 | cut -f1 -d"."`

if [ ${S1_hap} == 1 ];
then
    S1_HPRC_hap="hap2"
else
    S1_HPRC_hap="hap1"
fi

if [ ${S2_hap} == 1 ];
then
    S2_HPRC_hap="hap2"
else
    S2_HPRC_hap="hap1"
fi

FASTA1=${FASTADIR}${S1_ID}_${S1_HPRC_hap}/analysis/extract_hors_outputs/${S1_ID}_${S1_HPRC_hap}_${SAMPLE1}_chr12_hor_array.fasta

FASTA2=${FASTADIR}${S2_ID}_${S2_HPRC_hap}/analysis/extract_hors_outputs/${S2_ID}_${S2_HPRC_hap}_${SAMPLE2}_chr12_hor_array.fasta

echo "fasta 1:" $FASTA1
echo "fasta 2:" $FASTA2
TEMP_FASTA=${WORKDIR}/${SAMPLE1}_${SAMPLE2}.fa
cat $FASTA1 $FASTA2 > $TEMP_FASTA
time /private/home/mmastora/progs/centrolign/build/centrolign -v 3 --skip-calibration $TEMP_FASTA > $OUTDIR/pairwise_cigar_${SAMPLE1}_${SAMPLE2}.txt
rm $TEMP_FASTA
```
Trio only samples:
```
/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/fasta_list.all_sample_ids.in_nwk.trio_only.txt
```

Get pairwise distance matrix
```
for f in 100kb 50kb 500kb; do
  python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/cigar_to_distance.py /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps_${f}_flanks/chr12/all_pairs_pairwise/pairwise_cigar
done

for f in 100kb 50kb 500kb; do
  mv /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps_${f}_flanks/chr12/all_pairs_pairwise/pairwise_cigarpairwise_distance.csv /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps_${f}_flanks/chr12/all_pairs_pairwise/initial_test_nogaps_${f}_flanks_pairwise_distance.csv
done

for f in 100kb 50kb 500kb; do
  cp /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps_${f}_flanks/chr12/all_pairs_pairwise/initial_test_nogaps_${f}_flanks_pairwise_distance.csv /private/groups/patenlab/mira/
done
```

Combine distances matrices and infer tree
```
for f in 100kb 50kb 500kb; do
    python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/combine_HOR_flank_dist.py \
        -c /Users/miramastoras/Desktop/centrolign_all_pairs_flanks/initial_test_nogaps_${f}_flanks_pairwise_distance.csv \
        -f /Users/miramastoras/Desktop/distance_compare/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_hap.hD.txt \
        -o /Users/miramastoras/Desktop/centrolign_all_pairs_flanks/initial_test_nogaps_${f}_flanks_all_pairs_combined_dist
  done
```
Plot pairwise heatmap
```
conda activate tree_python

# combined distances as a spot check
for f in 100kb 50kb 500kb; do
  python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/centrolign_all_pairs_flanks/initial_test_nogaps_${f}_flanks_all_pairs_combined_dist_HOR_flank_dist_weighted.nwk \
        -s /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.txt \
        -p /Users/miramastoras/Desktop/centrolign_all_pairs_flanks/initial_test_nogaps_${f}_flanks_all_pairs_combined_dist_HOR_flank_dist_weighted.txt \
        -o /Users/miramastoras/Desktop/centrolign_all_pairs_flanks/combine_HOR_flank_dist_tree_and_aligment_dists_chr12_all_pairs_${f}_flanks
done

# just centrolign distances
for f in 100kb 50kb 500kb; do
  python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/centrolign_all_pairs_flanks/initial_test_nogaps_${f}_flanks_all_pairs_combined_dist_HOR_flank_dist_weighted.nwk \
        -s /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.txt \
        -p /Users/miramastoras/Desktop/centrolign_all_pairs_flanks/initial_test_nogaps_${f}_flanks_pairwise_distance.csv \
        -o /Users/miramastoras/Desktop/centrolign_all_pairs_flanks/combine_HOR_flank_dist_tree_all_pairs_alignment_dists_chr12_${f}_flanks
done
```
Comparing the trees
```
install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
library(phytools)
library(ape)

# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/chr12_HOR_flank_dist_weighted.nwk")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/centrolign_all_pairs_flanks/initial_test_nogaps_100kb_flanks_all_pairs_combined_dist_HOR_flank_dist_weighted.nwk")

samples <- readLines("/Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.txt")

obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/centrolign_all_pairs_flanks/HOR_vs_100kb_tree.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()
```
HOR vs 50kb tree
```
install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
library(phytools)
library(ape)

# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/chr12_HOR_flank_dist_weighted.nwk")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/centrolign_all_pairs_flanks/initial_test_nogaps_50kb_flanks_all_pairs_combined_dist_HOR_flank_dist_weighted.nwk")

samples <- readLines("/Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.txt")

obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/centrolign_all_pairs_flanks/HOR_vs_50kb_tree.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()
```
HOR vs 500kb tree
```
install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
library(phytools)
library(ape)

# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/chr12_HOR_flank_dist_weighted.nwk")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/centrolign_all_pairs_flanks/initial_test_nogaps_500kb_flanks_all_pairs_combined_dist_HOR_flank_dist_weighted.nwk")

samples <- readLines("/Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.txt")

obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/centrolign_all_pairs_flanks/HOR_vs_500kb_tree.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()
```
### Infer tree from centrolign 50kb, all pairs + 100 and then + 500kb of the flanks

```
for f in 100kb 500kb ; do
    python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/data_exploration/infer_tree.py /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps_${f}_flanks/chr12/all_pairs_pairwise/pairwise_cigar/ > /private/groups/patenlab/mira/centrolign/guide_tree_testing/method3_all_vs_all/initial_test_nogaps_${f}_flanks_chr12.all_pairs.nj.nwk
done

python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/data_exploration/infer_tree.py /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/all_pairs_pairwise/pairwise_cigar/ > /private/groups/patenlab/mira/centrolign/guide_tree_testing/method3_all_vs_all/initial_test_nogaps_chr12.all_pairs.nj.nwk
```

Plot it against centrolign alignment distances
```
for f in 100kb 500kb; do
  python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/all_pairs_flanks_tree/initial_test_nogaps_${f}_flanks_chr12.all_pairs.nj.nwk \
        -s /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.txt \
        -p /Users/miramastoras/Desktop/centrolign_all_pairs_flanks/initial_test_nogaps_${f}_flanks_pairwise_distance.csv \
        -o /Users/miramastoras/Desktop/all_pairs_flanks_tree/initial_test_nogaps_${f}_flanks_chr12.all_pairs.nj_vs_all_pairs_dists
done
```
compare with Sasha's latest HPRC trees

swap sample ids back to hprc convention
```
sed 's/\.1:/MIRA/g' /Users/miramastoras/Desktop/all_pairs_flanks_tree/initial_test_nogaps_500kb_flanks_chr12.all_pairs.nj.nwk | sed 's/\.2:/.1:/g' | sed 's/MIRA/.2:/g' > /Users/miramastoras/Desktop/all_pairs_flanks_tree/initial_test_nogaps_500kb_flanks_chr12.all_pairs.nj.HPRC_naming.nwk

sed 's/\.1:/MIRA/g' /Users/miramastoras/Desktop/all_pairs_flanks_tree/initial_test_nogaps_100kb_flanks_chr12.all_pairs.nj.nwk | sed 's/\.2:/.1:/g' | sed 's/MIRA/.2:/g' > /Users/miramastoras/Desktop/all_pairs_flanks_tree/initial_test_nogaps_100kb_flanks_chr12.all_pairs.nj.HPRC_naming.nwk

sed 's/\.1:/MIRA/g' /Users/miramastoras/Desktop/centrolign_all_pairs/chr12_centrolign_all_pairs_nj_tree.nwk | sed 's/\.2:/.1:/g' | sed 's/MIRA/.2:/g' > /Users/miramastoras/Desktop/centrolign_all_pairs/chr12_centrolign_all_pairs.nj.HPRC_naming.nwk

sed 's/\.1/MIRA/g' /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.txt | sed 's/\.2/.1/g' | sed 's/MIRA/.2/g' > /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.HPRC_naming.txt
```

```
install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
library(phytools)
library(ape)

# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/Sasha_HPRC_trees_2_25_25/HPRC_chr12_P_Q_mira.2.25.25.upgma.nwk")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/all_pairs_flanks_tree/initial_test_nogaps_500kb_flanks_chr12.all_pairs.nj.HPRC_naming.nwk")

samples <- readLines("/Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.HPRC_naming.txt")

obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/all_pairs_flanks_tree/HPRC_chr12_P_Q_mira.2.25.25.upgma_vs_500kb_flanks_all_pairs_nj.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()
```
100kb tree
```
install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
library(phytools)
library(ape)

# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/Sasha_HPRC_trees_2_25_25/HPRC_chr12_P_Q_mira.2.25.25.upgma.nwk")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/all_pairs_flanks_tree/initial_test_nogaps_100kb_flanks_chr12.all_pairs.nj.HPRC_naming.nwk")

samples <- readLines("/Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.txt")

obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/all_pairs_flanks_tree/HPRC_chr12_P_Q_mira.2.25.25.upgma_vs_100kb_flanks_all_pairs_nj.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()
```
HOR tree
```
install.packages("phytools",repos="https://cloud.r-project.org",quiet=TRUE)
library(phytools)
library(ape)

# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/Sasha_HPRC_trees_2_25_25/HPRC_chr12_P_Q_mira.2.25.25.upgma.nwk")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/centrolign_all_pairs/chr12_centrolign_all_pairs.nj.HPRC_naming.nwk")

samples <- readLines("/Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.txt")

obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/all_pairs_flanks_tree/HPRC_chr12_P_Q_mira.2.25.25.upgma_vs_all_pairs_nj.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()
```

Sasha's nj trees
```
# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/Sasha_HPRC_trees_2_25_25/HPRC_chr12_P_Q_mira.2.25.25.rnj.format5.nwk")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/all_pairs_flanks_tree/initial_test_nogaps_100kb_flanks_chr12.all_pairs.nj.HPRC_naming.nwk")

samples <- readLines("/Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.txt")

obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/all_pairs_flanks_tree/HPRC_chr12_P_Q_mira.2.25.25.rnj_vs_100kb_flanks_all_pairs_nj.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()
```
Tree using 100kb flanks only vs tree using upgma distances and centrolign distances in the formula
```
# create trees
original_tree <- ape::read.tree("/Users/miramastoras/Desktop/all_pairs_flanks_tree/initial_test_nogaps_100kb_flanks_chr12.all_pairs.nj.nwk")
new_tree <- ape::read.tree("/Users/miramastoras/Desktop/centrolign_all_pairs_flanks/initial_test_nogaps_100kb_flanks_all_pairs_combined_dist_HOR_flank_dist_weighted.nwk")

samples <- readLines("/Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.txt")

obj<-cophylo(original_tree,new_tree)
svg(filename="/Users/miramastoras/Desktop/all_pairs_flanks_tree/centrolign_100kb_flanks_all_pairs_only_vs_weighted_sum_old_tree.svg")
plot(obj, pts=FALSE,fsize=c(0.15,0.15),link.type="curved",link.lty="solid",link.col=make.transparent("blue",0.25),part=0.44)
dev.off()
```
