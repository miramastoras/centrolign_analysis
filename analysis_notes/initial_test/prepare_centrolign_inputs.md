### Preparing inputs to centrolign for initial test on chr12

Using 190 hprc release 2 assemblies, annotated by Julian

#### 1. Align each haplotype to CHM13

wdl: https://github.com/miramastoras/hpp_production_workflows/blob/master/QC/wdl/tasks/Asm2AsmAlignerPaf.wdl

batch submission: https://docs.google.com/spreadsheets/d/1is_jiWsDoqj_1QIcunGJLoojmToX3z9TFvTSEWCY2jA/edit?gid=0#gid=0

#### 2. Extract chr12 HOR arrays for each assembly and format correctly

Jordan's script [locate_hors_from_censat.py](https://github.com/jeizenga/centromere-scripts/blob/main/benchmarking/locate_hors_from_censat.py) takes in alignment to chm13, censat annotation and creates a single bed entry for every hor array, filtering out hor arrays that are mapped to the wrong chromosome, or any that are broken / discontiguous. It also only keeps cross-chromosomal groups where all of them are fully resolved because in these cases if we have loose contigs that we can’t assign to a chromosome, we don’t know for sure that they DIDN’T come from one of the ones that look complete. This could potentially result in low sample sizes for cross-chromosomal groups containing both acrocentrics and metacentrics  

Prepare inputs to script:
https://docs.google.com/spreadsheets/d/1is_jiWsDoqj_1QIcunGJLoojmToX3z9TFvTSEWCY2jA/edit?gid=1739270624#gid=1739270624

list censat files for csv
```
ls | while read line ; do realpath ${line}/*.as_hor_sf.bed ; done > /private/groups/patenlab/mira/centrolign/annotations/as_hor_files.csv

ls | while read line ; do realpath ${line}/*.cenSat.bed ; done > /private/groups/patenlab/mira/centrolign/annotations/cenSat_files.csv
```

Run sbatch script

https://github.com/miramastoras/centrolign_analysis/tree/main/batch_submissions/extract_hors_from_assemblies_sbatch

#### 3. Run centrolign

List chr12 fasta files
```
ls | while read line ; do realpath $line/*.fasta ; done | grep "hor_array" > fasta_list.txt

# 161 / 189 did not have chr12 filtered out
```
Combine files
```
cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_from_assemblies_sbatch

cat fasta_list.txt | while read line ; do cat $line ; done > chr12_hprc_r2_initial_test.fasta
samtools faidx chr12_hprc_r2_initial_test.fasta
```
Remove samples which aren't in the tree
```
# get list of sample names as they'd be listed in tree
cat /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_from_assemblies_sbatch/fasta_list.txt | while read line ; do basename $line | cut -f3 -d"_" ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/chr12_initial_test_all_samples.txt

# print those in tree
SAMPLES=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/chr12_initial_test_all_samples.txt
NWK=/private/groups/patenlab/jeizenga/centromere/chr12/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.nwk.txt

while IFS= read -r pattern; do
  if grep -q "$pattern" $NWK; then
    echo "$pattern"
  fi
done < $SAMPLES | wc -l

# 128 / 161

while IFS= read -r pattern; do
  if grep -q "$pattern" $NWK; then
    echo "$pattern"
  fi
done < $SAMPLES > chr12_hprc_r2_initial_test_in_nwk.txt
```

Make new fasta containing only samples found in tree
```
cat chr12_hprc_r2_initial_test_in_nwk.txt | while read line ; do
    grep $line /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_from_assemblies_sbatch/fasta_list.txt
  done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/fasta_list_inside_nwk.txt

# combine fastas
cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/fasta_list_inside_nwk.txt | while read line ; do cat $line ; done > chr12_hprc_r2_initial_test_inside_tree.fasta
samtools faidx chr12_hprc_r2_initial_test_inside_tree.fasta
```

Run centrolign
```
#!/bin/bash
#SBATCH --job-name=centrolign_initial_chr12
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=800gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign_%A.log
#SBATCH --time=72:00:00

time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/jobstore/ \
    -T /private/groups/patenlab/jeizenga/centromere/chr12/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.nwk.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/chr12_hprc_r2_initial_test_inside_tree.fasta \
    > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/chr12_hprc_r2_initial_test_inside_tree.centrolign.gfa
```

Do another run, with subset of 50 samples
```
cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/chr12_hprc_r2_initial_test_in_nwk.txt | head -n 50 | while read line ; do
    grep $line /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_from_assemblies_sbatch/fasta_list.txt
  done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_50/fasta_list_inside_nwk.txt

# combine fastas
cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_50/fasta_list_inside_nwk.txt | while read line ; do cat $line ; done > chr12_hprc_r2_initial_test_inside_tree.first50.fasta
samtools faidx chr12_hprc_r2_initial_test_inside_tree.first50.fasta
```

```
#!/bin/bash
#SBATCH --job-name=centrolign_initial_chr12_first50
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=800gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign_%A.log
#SBATCH --time=72:00:00

time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_50/jobstore/ \
    -T /private/groups/patenlab/jeizenga/centromere/chr12/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.nwk.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_50/chr12_hprc_r2_initial_test_inside_tree.first50.fasta \
    > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_50/chr12_hprc_r2_initial_test_inside_tree.first50.centrolign.gfa
```

### Preparing inputs to centrolign for testing on chrX,Y,10,6,17

#### 1. Use new wdl to run pre-processing procedure for all chromosomes

batch submission:
https://github.com/miramastoras/centrolign_analysis/tree/main/batch_submissions/extract_hors/initial_test

#### 2. Prepare centrolign inputs


For chr6,10,17 need to remove cenhap label from tree
```
sed 's/_ch[1-9][0-9]*//g' chr6_HPRC_cenhap.nwk > chr6_HPRC.nwk
sed 's/_ch[1-9][0-9]*//g' chr17_HPRC_cenhap.nwk > chr17_HPRC.nwk
sed 's/_ch[1-9][0-9]*//g' chr10_HPRC_cenhap.nwk > chr10_HPRC.nwk
```

Get lists of fasta files containing complete HOR for each chrom
```
cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors/initial_test

for CHR in chrY chrX chr10 chr6 chr17 ; do
    mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/${CHR}/
    ls | grep hap | \
    while read line ; do
        realpath $line/analysis/extract_hors_outputs/*.fasta
      done | grep ${CHR} > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/${CHR}/fasta_list.txt
  done
```

Get list of samples that are in the nwk tree

https://docs.google.com/spreadsheets/d/1is_jiWsDoqj_1QIcunGJLoojmToX3z9TFvTSEWCY2jA/edit?gid=452898455#gid=452898455

we don't have a tree for chrY.

```
# get list of sample names as they'd be listed in tree
for CHR in chrY chrX chr10 chr6 chr17 ; do  
    cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/${CHR}/fasta_list.txt | while read line ; do basename $line | cut -f2-3 -d"." ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/${CHR}/fasta_list.all_sample_ids.txt
  done

# chrY chrX
for CHR in chr10 chr6 chr17 ; do   
    SAMPLES=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/${CHR}/fasta_list.all_sample_ids.txt
    NWK=`grep ${CHR} /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/batch_submissions/centrolign/initial_test_sbatch/centrolign_initial_test_sbatch.csv | cut -f 3 -d","`
    while IFS= read -r pattern; do
      if grep -q "$pattern" $NWK; then
        echo "$pattern
        fi
    done < $SAMPLES > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/${CHR}/fasta_list.all_sample_ids.in_nwk.txt
  done

# combine fastas
for CHR in chr10 chr6 chr17 ; do  
    cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/${CHR}/fasta_list.all_sample_ids.in_nwk.txt | while read line ; do
    grep $line /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/${CHR}/fasta_list.txt
    done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/${CHR}/fasta_list.inside_nwk.txt
done

# combine fastas
# for CHR in chrY chrX  ; do  
for CHR in chr10 chr6 chr17 ; do
    cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/${CHR}/fasta_list.inside_nwk.txt | while read line ; do cat $line ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/${CHR}/initial_test_${CHR}.fasta
    samtools faidx /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/${CHR}/initial_test_${CHR}.fasta
done
```
List for csv: https://docs.google.com/spreadsheets/d/1is_jiWsDoqj_1QIcunGJLoojmToX3z9TFvTSEWCY2jA/edit?gid=452898455#gid=452898455

```
for CHR in chr6 chr10 chr17 ; do  
    realpath /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/${CHR}/initial_test_${CHR}.fasta
  done
```
Run as a batch submission:


Deal with chrX separately:

chrX newick file is labelled as all males but contains female samples, and all samples have .0 instead of .1 or .2. For now just going to run the male samples  
```
# remove haplotype label from samples containing chrX
cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/chrX/fasta_list.txt | while read line ; do basename $line | cut -f2 -d"." | sort | uniq ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/chrX/fasta_list.all_sample_ids.txt

# grep sample names that are in the guide tree
for CHR in chrX ; do   
    SAMPLES=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/${CHR}/fasta_list.all_sample_ids.txt
    NWK=/private/groups/patenlab/jeizenga/centromere/chrX/KGP4_Hg38_MAC10_chrX_MALES_chm13_HG002_HG005_HuRef.nwk.txt
    while IFS= read -r pattern; do
      if grep -q "$pattern" $NWK; then
        echo "$pattern"
        fi
    done < $SAMPLES > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/${CHR}/fasta_list.all_sample_ids.in_nwk.txt
done

# remove female samples
grep "female" /private/groups/patenlab/mira/centrolign/annotations/guide_trees/HPRC_Intermediate_Assembly_Data_Status.csv | cut -f 2 -d"," > /private/groups/patenlab/mira/centrolign/annotations/guide_trees/HPRC_Intermediate_Assembly_Data_Status.female_sample_ids.txt

grep -Fvxf /private/groups/patenlab/mira/centrolign/annotations/guide_trees/HPRC_Intermediate_Assembly_Data_Status.female_sample_ids.txt /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/chrX/fasta_list.all_sample_ids.in_nwk.txt > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/chrX/fasta_list.all_sample_ids.in_nwk.male_only.txt

# replace naming inside fasta file to match newick
cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/chrX/fasta_list.all_sample_ids.in_nwk.male_only.txt | while read line ; do
    sed 's/${line}.2//g'


```
