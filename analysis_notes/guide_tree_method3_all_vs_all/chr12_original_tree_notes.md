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

python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/chr12_HOR_flank_dist_weighted.nwk \
        -s /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.txt \
        -p /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/chr12_HOR_flank_dist_weighted.txt \
        -o /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/combine_HOR_flank_dist_chr12_all_pairs_

#
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap.py \
        -t /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/chr12_HOR_flank_dist_weighted.nwk \
        -s /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/fasta_list.all_sample_ids.in_nwk.trio_only.txt \
        -p /Users/miramastoras/Desktop/centrolign_all_pairs/pairwise_distance.csv \
        -o /Users/miramastoras/Desktop/all_pairs_weighted_sum_trio_ch12/combine_HOR_flank_dist_chr12_all_pairs_aln
```

### Test expanding out sequence into the flanks for chr12

#### 100 kb flanks

#### 50 kb flanks

#### 500 kb flanks 
