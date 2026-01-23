- for every pair:
- Take CDR coordinates for each sample
- calculate mutation rate in 1,000 bp windows inside CDR
- randomly shuffle 10,000 windows? across rest of the array, calculate mutation rate
- compare the distribution 


mkdir permuted_beds


for i in {1..10000}; do
  bedtools shuffle \
    -i cdr_windows.bed \
    -g genome.sizes \
    -incl centromere_minus_CDR.bed \
    -noOverlapping \
    > permuted_beds/perm_${i}.bed
done
