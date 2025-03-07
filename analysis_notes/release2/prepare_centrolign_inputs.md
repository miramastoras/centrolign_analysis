# Preparing inputs to centrolign

Aligned to CHM13, then extracted HOR arrays that are contiguous
https://docs.google.com/spreadsheets/d/1is_jiWsDoqj_1QIcunGJLoojmToX3z9TFvTSEWCY2jA/edit?gid=1520027858#gid=1520027858

https://github.com/miramastoras/centrolign_analysis/tree/main/batch_submissions/Asm2AsmAlignerPaf/release2
https://github.com/miramastoras/centrolign_analysis/tree/main/batch_submissions/extract_hors_HPRC/release2

For each chr, get list of samples with contiguous arrays
```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/contiguous_HORs/

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

for chr in "${chromosomes[@]}"
do
    echo "Processing $chr"
    ls | grep hap | \
    while read line ; do
        ls $line/analysis/extract_hors_HPRC_outputs/*${chr}_hor_array.fasta | cut -f4 -d"/" | cut -f6 -d"_"
      done | grep -v "cannot access" > /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/contiguous_HORs/HPRC_release2_contiguous_HORs_${chr}.txt
done

# Count number of contiguous arrays per chromosome
for chr in "${chromosomes[@]}"
do  
  count=`wc -l /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/contiguous_HORs/HPRC_release2_contiguous_HORs_${chr}.txt | cut -f1 -d" "`
  echo $chr $count
done
```
Parse logs for reason for filtering
```
cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/slurm_logs

for log in *.log ;
  do grep "filtering" $log | grep -v "\["
done > filter_log.txt

# 491 filtering for gaps
# 2013 	filtering for discontiguous array
# 7 filtering for discordant mapped chromosome chr18
# 199 filtering HOR group for incomplete assignments
# 809 filtering for discontiguous array over chromosome _ in a resolved HOR group
# 10 filtering for being at end of contig of length _ in unique type
# 14 filtering for being at end of contig of length _ in a resolved HOR group

grep -v "filtering for discontiguous array over chromosome" filter_log.txt | grep -v "filtering for being at end of contig of length" | grep -v "filtering HOR group for incomplete assignments" | grep -v "filtering for discordant mapped chromosome" | grep -v "filtering for discontiguous array" | grep -v "filtering for gaps"

for log in *.log ;
  do grep "filtering for discordant HOR type" $log | grep -v "\["
done > filter_log2.txt
```
Get list of fasta paths per chromosome
```sh

cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

for chr in "${chromosomes[@]}"
do
    echo "Processing $chr"
    ls | grep hap | \
    while read line ; do
        echo $line/analysis/extract_hors_HPRC_outputs/*.fasta | grep ${chr}
      done
done
```
