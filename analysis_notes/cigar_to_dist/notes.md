## Comparing different formulas for converting a cigar string to distance

The purpose of this analysis is to find the best way to convert the cigar string to a pairwise distance for centrolign.

Formulas tested:

1. Jordan's original formula

```
1.0 - (2.0 * matches)
       ---------------
     (ref_len + query_len)
         |          |_ matches + mismatches + unaligned bases in query
         |_ matches + mismatches + unaligned bases in ref
```
2. Jordan's proposal involving the proportion of unaligned sequence:
```
identity:
(proportion aligned) * (matches / (matches + mismatches))

distance:
(proportion unaligned) * (mismatches / matches + mismatches)
```
3. Benedict's proposal:
```
mismatches/(matches+mismatches)
```
> We decided to scrap this one, because it breaks down in cases like this `210=2359579D3069392I168=` where the two sequences have a very small # of aligned bases

### 1. Implementing the formulas

Selecting two test cigar strings to operate on:

```
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/cigar_to_distance.py \
  -a /Users/miramastoras/Desktop/github_repos/centrolign_analysis/analysis_notes/cigar_to_dist/pairwise_cigar1/ \
  -d 1
```

```
Cigar: 210=2359579D3069392I168=

Proportion aligned = (ref aligned + query aligned / ref total + query total)
1 - ((proportion aligned) * (matches / (mismatches + matches)))
1 - ( (  (378*2) / ((378*2)+2359579+3069392+ ) * (378 / (0 + 378))) = 0.999860766

Cigar: 915874=6X100=
Proportion aligned = (ref aligned + query aligned / ref total + query total)
Proportion aligned = 1
1 - ((proportion aligned) * (matches / (mismatches + matches)))
1 - (1 * (915974 / 915980 ) ) = 0.00000655036

Cigar: 100=2359D3060I1689=
Proportion aligned = (ref aligned + query aligned / ref total + query total)
Proportion aligned = ((100+1689)*2) / (((100+1689)*2) + 2359+3060) = 0.39768811826

1 - ((proportion aligned) * (matches / (mismatches + matches)))
1 - ((0.39768811826) * (100+1689/100+1689) = 0.60231188174

Cigar: 100=2359I16890=
Proportion aligned = (ref aligned + query aligned / ref total + query total)
Proportion aligned = ((16890+100) *2) / (((16890+100) *2) + 2359) = 0.93508351908

1 - ((proportion aligned) * (matches / (mismatches + matches)))
1 - ((0.93508351908) * (16990 / 16990)) = 0.06491648092
```
