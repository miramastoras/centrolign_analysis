### Implementing Jordan's idea to assess alignability

- read in list of samples
- for each sample, get length of asat array
- Created mapping of smush coordinate system from 0 to 1 to original coordinate system. dict key=original, value = smush
- load in cigar strings from each sample
- for every sample:
- loop through every cigar string containing that sample. if sample name is listed second in file name, rever
- For all cigar strings from that sample: take all cigar strings with this sample as reference. At every base in this sample, tabulate the proportion of bases aligned (mismatch or match) or unaligned (I/D)
-  draw 100 curves on 1 plot, can also get a mean curve
