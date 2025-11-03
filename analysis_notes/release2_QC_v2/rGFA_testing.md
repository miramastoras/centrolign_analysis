## Create simple centrolign graph with flanks for Glenn to test rGFA decomposition

Using test set of 50 samples with flanks
```sh
# get all pairs distances
while read -r s1; do
    while read -r s2; do
        [[ "$s1" < "$s2" ]] && grep -w $s1 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/chr12_r2_QC_v2_centrolign_pairwise_distance.csv | grep -w $s2 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/chr12_r2_QC_v2_centrolign_pairwise_distance.csv
    done < /private/groups/patenlab/mira/centrolign/batch_submissions/runtime_debugging/r2_QC_v2_flanks_test/chr12/sample_names.txt
done < /private/groups/patenlab/mira/centrolign/batch_submissions/runtime_debugging/r2_QC_v2_flanks_test/chr12/sample_names.txt | sort -n -t ',' -k 3 | uniq > /private/groups/patenlab/mira/centrolign/batch_submissions/runtime_debugging/r2_QC_v2_flanks_test/chr12/pairwise_distances.txt


ls | grep fasta | while read line ; do cut -f1 -d"_" | while read line ; do grep -w $line /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/chr12_r2_QC_v2_centrolign_pairwise_distance.csv | grep -w "HG04228.2" ; done ;done


# select 3 relatively closely related samples
HG01175.1,HG04228.2,NA19682.1

# pairwise dists
HG01175.1,HG04228.2,0.38270151701245880282
HG01175.1,NA19682.1,0.58092091692202574293
HG04228.2,NA19682.1,0.53691252325430671721
```
Combine fastas
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/rGFA_tests/chr12_3_samples_100kb_flank/

cat /private/groups/patenlab/mira/centrolign/batch_submissions/runtime_debugging/r2_QC_v2_flanks_test/chr12/extract_fastas_100kb_flank/HG01175.1_chr12_hor_array.fasta /private/groups/patenlab/mira/centrolign/batch_submissions/runtime_debugging/r2_QC_v2_flanks_test/chr12/extract_fastas_100kb_flank/HG04228.2_chr12_hor_array.fasta /private/groups/patenlab/mira/centrolign/batch_submissions/runtime_debugging/r2_QC_v2_flanks_test/chr12/extract_fastas_100kb_flank/NA19682.1_chr12_hor_array.fasta > /private/groups/patenlab/mira/centrolign/rGFA_tests/chr12_3_samples_100kb_flank/chr12_3smps_100kb_flank.fasta


# set up directory
mkdir -p /private/groups/patenlab/mira/centrolign/rGFA_tests/chr12_3_samples_100kb_flank/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/rGFA_tests/chr12_3_samples_100kb_flank/logs
```

```sh
#!/bin/bash
#SBATCH --job-name=chr12_3_samples_rGFA_test
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=400gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=12:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/rGFA_tests/chr12_3_samples_100kb_flank/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr12_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  /private/groups/patenlab/mira/centrolign/rGFA_tests/chr12_3_samples_100kb_flank/chr12_3smps_100kb_flank.fasta \
  > /private/groups/patenlab/mira/centrolign/rGFA_tests/chr12_3_samples_100kb_flank/chr12.3smps_100kb_flank.centrolign.gfa
```
