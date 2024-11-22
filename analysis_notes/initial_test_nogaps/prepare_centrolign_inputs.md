### Check how many hor arrays had gaps per sample

```
cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors/initial_test_nogaps/slurm_logs

ls | while read line ; do
    sample=`grep "SAMPLE_ID=" $line | cut -f2 -d"="`
    gaps=`grep "filtering for gaps:" $line | grep -v "\[" | wc -l`
    echo $sample,$gaps
done > gaps_per_sample.txt

# count of hor arrays filtered
awk -F',' '{sum+=$2;} END{print sum;}' gaps_per_sample.txt  
# 195

# distribution per sample
cut -f2 -d"," gaps_per_sample.txt | sort | uniq -c

68 0
69 1
36 2
14 3
3 4
```

## prepare inputs to centrolign

Get lists of fasta files containing complete HOR for each chrom
```
cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors/initial_test_nogaps

for CHR in chr10 chr6 chr17 ; do
    mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/${CHR}/
    ls | grep hap | \
    while read line ; do
        realpath $line/analysis/extract_hors_outputs/*.fasta
      done | grep ${CHR} > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/${CHR}/fasta_list.txt
  done
```

Get list of samples that are in the nwk tree

https://docs.google.com/spreadsheets/d/1is_jiWsDoqj_1QIcunGJLoojmToX3z9TFvTSEWCY2jA/edit?gid=452898455#gid=452898455

we don't have a tree for chrY.

```
# get list of sample names as they'd be listed in tree
for CHR in chr10 chr6 chr17 ; do  
    cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/${CHR}/fasta_list.txt | while read line ; do basename $line | cut -f3 -d"_" ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/${CHR}/fasta_list.all_sample_ids.txt
  done

# chrY chrX
for CHR in chr10 chr6 chr17 ; do   
    SAMPLES=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/${CHR}/fasta_list.all_sample_ids.txt
    NWK=`grep ${CHR} /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/batch_submissions/centrolign/initial_test_sbatch/centrolign_initial_test_sbatch.csv | cut -f 3 -d","`
    while IFS= read -r pattern; do
      if grep -q "$pattern" $NWK; then
        echo "$pattern"
        fi
    done < $SAMPLES > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/${CHR}/fasta_list.all_sample_ids.in_nwk.txt
  done

# combine fastas
for CHR in chr10 chr6 chr17 ; do  
    cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/${CHR}/fasta_list.all_sample_ids.in_nwk.txt | while read line ; do
    grep $line /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/${CHR}/fasta_list.txt
    done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/${CHR}/fasta_list.inside_nwk.txt
done

# combine fastas
# for CHR in chrY chrX  ; do  
for CHR in chr10 chr6 chr17 ; do
    cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/${CHR}/fasta_list.inside_nwk.txt | while read line ; do cat $line ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/${CHR}/initial_test_${CHR}.fasta
    samtools faidx /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/${CHR}/initial_test_${CHR}.fasta
done
```
List for csv: https://docs.google.com/spreadsheets/d/1is_jiWsDoqj_1QIcunGJLoojmToX3z9TFvTSEWCY2jA/edit?gid=452898455#gid=452898455

```
for CHR in chr6 chr10 chr17 ; do  
    realpath /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/${CHR}/initial_test_${CHR}.fasta
  done
```

Chr X
```
cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chrX/fasta_list.txt | while read line ; do basename $line | cut -f3 -d"_" | sort | uniq ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chrX/fasta_list.all_sample_ids.txt

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
    sed "s|${line}.1|${line}.0|g" /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors/initial_test/${line}_hap2/analysis/extract_hors_outputs/${line}_hap2.${line}.1.chrX_hor_array.fasta > /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors/initial_test/${line}_hap2/analysis/extract_hors_outputs/${line}_hap2.${line}.0.chrX_hor_array.fasta
  done

# get fasta locations
cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/chrX/fasta_list.all_sample_ids.in_nwk.male_only.txt | while read line ; do
ls /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors/initial_test/${line}_hap2/analysis/extract_hors_outputs/${line}_hap2.${line}.0.chrX_hor_array.fasta
done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/chrX/all_sample_ids.in_nwk.male_only.full_path.txt

# combine fastas
cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/chrX/all_sample_ids.in_nwk.male_only.full_path.txt | while read line ; do cat $line ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/chrX/initial_test_chrX.male_only.fasta
samtools faidx /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/chrX/initial_test_chrX.male_only.fasta
```
