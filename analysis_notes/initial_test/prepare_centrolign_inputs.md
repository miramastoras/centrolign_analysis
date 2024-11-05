### Preparing inputs to centrolign for initial test

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
