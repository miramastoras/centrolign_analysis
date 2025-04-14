## Centrolign analysis of chr 6 for release 2

### 0. Prepare data

#### Cenhap distances

Cenhap guide trees and distance matrix provided by Chuck and Sasha Langley on 4/11/2025. These 7 samples are removed due to large deletions:
```
NA20827.1
HG00642.2
HG00344.2
HG00639.1
HG02622.1
HG02145.1
NA20799.2
```
distance matrix from p and q
```sh
/private/groups/patenlab/mira/centrolign/annotations/guide_trees/HPRC_chr6_58200000_58286706_61058390_61121735_het68_m_mira_dgp_rnj.m
```
Reformat matrix file to enable reading into python script
```
# convert spaces to comma, remove trailing commas
sed 's/ /,/g' HPRC_chr6_58200000_58286706_61058390_61121735_het68_m_mira_dgp_rnj.m
 | sed 's/,$//' | sed 's/\t//' | tail -n +2 > matrix_no_header.csv

# convert row names to column names
N=384
cut -f1 -d"," matrix_no_header.csv | sed 's/$/,/g' | tr -d '\n' | cut -f1-$N -d"," | sed 's/^/0,/g' > rowname_to_colname.csv

cat rowname_to_colname.csv matrix_no_header.csv > HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj.reformatted.m

rm matrix_no_header.csv rowname_to_colname.csv
```
