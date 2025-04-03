## Centrolign analysis of chr 12 for release 2

- Using chr12 as a first pass to implement our analysis plan, before automating it for the other chroms.

### 0. Prepare data

#### Cenhap distances

Cenhap guide trees and distance matrix provided by Chuck and Sasha Langley on 4/2/2025, containing all the samples with contiguous HORs from release 2, except HG00741.1 because it has a large deletion(?) relative to chm13 on the q side of the array

```sh
/private/groups/patenlab/mira/centrolign/annotations/guide_trees/HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj.m
/private/groups/patenlab/mira/centrolign/annotations/guide_trees/HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj.nwk
/private/groups/patenlab/mira/centrolign/annotations/guide_trees/HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj_upgma.nwk
```

Reformat matrix file to enable reading into python script
```
# convert spaces to comma, remove trailing commas
sed 's/ /,/g' HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj.m | sed 's/,$//' | sed 's/\t//' | tail -n +2 > matrix_no_header.csv

# convert row names to column names
N=384
cut -f1 -d"," matrix_no_header.csv | sed 's/$/,/g' | tr -d '\n' | cut -f1-$N -d"," | sed 's/^/0,/g' > rowname_to_colname.csv

cat rowname_to_colname.csv matrix_no_header.csv > HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj.reformatted.m

rm matrix_no_header.csv rowname_to_colname.csv
```
#### Centrolign all pairs distances

Convert cigar to distance matrix
```
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/cigar_to_distance.py /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr12/pairwise_cigar/
```

### 1. Create refined tree combining centrolign all pairs distances and cenhap distances

```
conda activate skbio

python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/combine_HOR_flank_dist.py \
    -c /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/chr12/pairwise_distance.csv \
    -f /private/groups/patenlab/mira/centrolign/annotations/guide_trees/HPRC_chr12_34543495_34593492_37202490_37286092_het66.7_m_mira_dgp_rnj.reformatted.m \
    -o /private/groups/patenlab/mira/centrolign/guide_tree_testing/release2_weighted_sum/HPRC_r2_chr12_cenhap_20250402_centrolign_all_pairs
```
