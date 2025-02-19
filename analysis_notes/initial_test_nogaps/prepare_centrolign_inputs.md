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

for CHR in chr10 chr6 chr17 chr12; do
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
for CHR in chr10 chr6 chr17 chr12 ; do  
    cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/${CHR}/fasta_list.txt | while read line ; do basename $line | cut -f3 -d"_" ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/${CHR}/fasta_list.all_sample_ids.txt
  done

# chrY chrX
for CHR in chr10 chr6 chr17 chr12 ; do   
    SAMPLES=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/${CHR}/fasta_list.all_sample_ids.txt
    NWK=`grep ${CHR} /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/batch_submissions/centrolign/initial_test_sbatch/centrolign_initial_test_sbatch.csv | cut -f 3 -d","`
    while IFS= read -r pattern; do
      if grep -q "$pattern" $NWK; then
        echo "$pattern"
        fi
    done < $SAMPLES > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/${CHR}/fasta_list.all_sample_ids.in_nwk.txt
  done

# combine fastas
for CHR in chr10 chr6 chr17 chr12; do  
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

for CHR in chr12; do
  cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/${CHR}/fasta_list.inside_nwk.txt | while read line ; do cat $line ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/${CHR}/initial_test_no_gaps_${CHR}.fasta
  samtools faidx /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/${CHR}/initial_test_no_gaps_${CHR}.fasta
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

Run centrolign, chr12, outputting pairwise cigars.
```
#!/bin/bash
#SBATCH --job-name=centrolign_chr12_nogaps_all_samples_pairwise
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=500gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign.log
#SBATCH --time=7-00:00


/private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/jobstore \
    -T /private/groups/patenlab/jeizenga/centromere/chr12/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.nwk.txt \
    -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/pairwise_cigars/ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/initial_test_no_gaps_chr12.fasta \
    > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/initial_test_no_gaps_chr12.centrolign.gfa
```
Submit centrolign for chr17 - submitting initial alignment first, then restart using multithreaded approach to output pairwise cigar
```
#!/bin/bash
#SBATCH --job-name=centrolign_chr17_nogaps
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=500gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign.log
#SBATCH --time=7-00:00

/private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr17/jobstore \
    -T /private/groups/patenlab/mira/centrolign/annotations/guide_trees/chr17_HPRC.nwk \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr17/initial_test_chr17.fasta \
    > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr17/initial_test_no_gaps_chr17.centrolign.gfa
```
```
-A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr17/pairwise_cigars/pairwise_cigar
```
chr 10
```
#!/bin/bash
#SBATCH --job-name=centrolign_chr10_nogaps
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=500gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign.log
#SBATCH --time=7-00:00

/private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr10/jobstore \
    -T /private/groups/patenlab/mira/centrolign/annotations/guide_trees/chr10_HPRC.nwk \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr10/initial_test_chr10.fasta \
    > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr10/initial_test_no_gaps_chr10.centrolign.gfa
```

chr6
```
#!/bin/bash
#SBATCH --job-name=centrolign_chr6_nogaps
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=500gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign.log
#SBATCH --time=7-00:00

/private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr6/jobstore \
    -T /private/groups/patenlab/mira/centrolign/annotations/guide_trees/chr6_HPRC.nwk \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr6/initial_test_chr6.fasta \
    > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr6/initial_test_no_gaps_chr6.centrolign.gfa
```
#### Testing newer chr12 guide tree provided by Sasha

```
#!/bin/bash
#SBATCH --job-name=centrolign_chr12_cenhap_update
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=500gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign.log
#SBATCH --time=7-00:00


/private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12_cenhap_update/jobstore/ \
    -T /private/groups/patenlab/mira/centrolign/annotations/guide_trees/chr12_cenhap_update.nwk \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/initial_test_no_gaps_chr12.fasta \
    -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/pairwise_cigars/ \
    > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12_cenhap_update/initial_test_no_gaps_chr12_cenhap_update.centrolign.gfa
```

### Potential sample haplotype swaps

confirmed by julian: HG01978, HG02257, HG03516

HG01784, HG02273, HG02451.2, NA19185.2

Check how many are trio samples:
```
grep "_mat_" /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors/initial_test_nogaps/extract_hors_initial_test_nogaps.csv >/private/groups/patenlab/mira/centrolign/test/samples_trio.txt
grep "_pat_" /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors/initial_test_nogaps/extract_hors_initial_test_nogaps.csv >> /private/groups/patenlab/mira/centrolign/test/samples_trio.txt

cut -f 1 /private/groups/patenlab/mira/centrolign/test/samples_trio.txt | cut -f 1 -d"_" | while read line ; do grep $line /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test_nogaps/chr12/fasta_list.all_sample_ids.in_nwk.txt ; done | sort | uniq | wc -l

# 107 / 125 samples are trio
```
