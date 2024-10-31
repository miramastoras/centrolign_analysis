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
ls | while read line ; do realpath $line/*.fasta ; done | grep "parnum" > fasta_list.txt

# 161 / 189 did not have chr12 filtered out
```
Combine files
```
cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_from_assemblies_sbatch

cat fasta_list.txt | while read line ; do cat $line ; done > chr12_hprc_r2_initial_test.fasta
```
Remove samples which aren't in the tree
```
/private/groups/patenlab/jeizenga/centromere/chr12/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.nwk.txt

cat fasta_list.txt | while read line ; do basename $line | cut -f1 -d"_" | sort | uniq ; done
```

Run centrolign
```
~/progs/centrolign/build/centrolign \
    -T /private/groups/patenlab/jeizenga/centromere/chr12/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.nwk.txt chr12_hprc_r2_initial_test.fasta > /private/groups/patenlab/mira/centrolign/initial_test/
```
