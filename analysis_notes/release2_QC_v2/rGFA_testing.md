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
### Create chr 12 graph with 150 samples, 100kb of flank, subsampled across clades


First check LDC I already have and color

```sh
docker run -u `id -u`:`id -g` \
    -v /private/groups:/private/groups/ \
    miramastoras/tree_heatmap:latest python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/pairwise_tree_heatmap_v2.py  \
    -t /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr12_r2_QC_v2_centrolign_all_pairs_nj_tree.format5.nwk \
   -s /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.chr12.samples.txt \
   -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/chr12_r2_QC_v2_centrolign_pairwise_distance.csv \
   -m "Centrolign all pairs distances" \
   -n "chr 12 NJ tree" \
   -d "All pairs Distances" \
   -o /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/plots/chr12_LDCs_ --no_labels \
   --highlight_samples /private/groups/patenlab/mira/centrolign/analysis/low_divergence_clades/max_dist_0.8_min_pairwise_0.95/chr12/_clades.csv
```

![LDCs for chr12](../plots/chr12_LDCs_pairwise_tree_heatmap.png)

Script to randomly shuffle 150 samples, giving equal samples for each clade, making sure to include CHM13.0, HG002.2 and HG002.1. Some samples don't have a clade, which is fine, just excluding those.

```py
#!/usr/bin/env python3
"""
Subsample to 150 samples equally across clades.
Always includes CHM13.0, HG002.1, HG002.2.
"""

import csv
import random
from collections import defaultdict

INPUT = "/private/groups/patenlab/mira/centrolign/analysis/low_divergence_clades/max_dist_0.8_min_pairwise_0.95/chr12/_clades.csv"
OUTPUT = "/private/groups/patenlab/mira/centrolign/analysis/low_divergence_clades/max_dist_0.8_min_pairwise_0.95/chr12/subsampled_150.txt"
TARGET = 150
REQUIRED = {"CHM13.0", "HG002.1", "HG002.2"}

random.seed(42)

# Load clade -> [samples]
clade_to_samples = defaultdict(list)
with open(INPUT) as f:
    reader = csv.DictReader(f)
    for row in reader:
        clade = row["Clade"].strip()
        sample = row["Sample"].strip()
        if clade == "Clade":  # skip header-like rows
            continue
        clade_to_samples[clade].append(sample)

# Start with required samples, track which clade they came from
selected = set(REQUIRED)

# Remove required samples from clade pools (avoid double-counting)
remaining = {}
for clade, samps in clade_to_samples.items():
    pool = [s for s in samps if s not in selected]
    random.shuffle(pool)
    remaining[clade] = pool

slots_left = TARGET - len(selected)

# Iteratively allocate equal slots across clades that still have samples
while slots_left > 0:
    active_clades = {c: pool for c, pool in remaining.items() if pool}
    if not active_clades:
        break

    per_clade = max(1, slots_left // len(active_clades))

    added_this_round = 0
    for clade in sorted(active_clades):
        take = min(per_clade, len(remaining[clade]), slots_left)
        if take <= 0:
            continue
        for s in remaining[clade][:take]:
            selected.add(s)
        remaining[clade] = remaining[clade][take:]
        slots_left -= take
        added_this_round += take
        if slots_left == 0:
            break

    if added_this_round == 0:
        break

# Build sample -> clade lookup
sample_to_clade = {}
for clade, samps in clade_to_samples.items():
    for s in samps:
        sample_to_clade[s] = clade

# Write output
with open(OUTPUT, "w") as f:
    for s in sorted(selected):
        clade = sample_to_clade.get(s, "unknown")
        f.write(f"{clade},{s}\n")

print(f"Selected {len(selected)} samples -> {OUTPUT}")

# Print per-clade breakdown
print("\nPer-clade counts in subsample:")
for clade in sorted(clade_to_samples):
    clade_samps = set(clade_to_samples[clade])
    n = len(selected & clade_samps)
    total = len(clade_samps)
    print(f"  {clade}: {n}/{total}")
```

Results: /private/groups/patenlab/mira/centrolign/analysis/low_divergence_clades/max_dist_0.8_min_pairwise_0.95/chr12/subsampled_150.txt

Plot just subsampled samples:

First check LDC I already have and color

```sh
docker run -u `id -u`:`id -g` \
    -v /private/groups:/private/groups/ \
    miramastoras/tree_heatmap:latest python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/pairwise_tree_heatmap_v2.py  \
    -t /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr12_r2_QC_v2_centrolign_all_pairs_nj_tree.format5.nwk \
   -s /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.chr12.samples.txt \
   -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/chr12_r2_QC_v2_centrolign_pairwise_distance.csv \
   -m "Centrolign all pairs distances" \
   -n "chr 12 NJ tree" \
   -d "All pairs Distances" \
   -o /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/plots/chr12_LDCs_subsampled_ --no_labels \
   --highlight_samples /private/groups/patenlab/mira/centrolign/analysis/low_divergence_clades/max_dist_0.8_min_pairwise_0.95/chr12/subsampled_150.txt
```
![LDCs for chr12](../plots/chr12_LDCs_subsampled_pairwise_tree_heatmap.png)


Now showing the plot with just the chr12 subsambled to 150:
```sh
docker run -u `id -u`:`id -g` \
    -v /private/groups:/private/groups/ \
    miramastoras/tree_heatmap:latest python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/pairwise_tree_heatmap_v2.py  \
    -t /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr12_r2_QC_v2_centrolign_all_pairs_nj_tree.format5.nwk \
   -s /private/groups/patenlab/mira/centrolign/analysis/low_divergence_clades/max_dist_0.8_min_pairwise_0.95/chr12/subsampled_150.sample_ids.txt \
   -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/chr12_r2_QC_v2_centrolign_pairwise_distance.csv \
   -m "Centrolign all pairs distances" \
   -n "chr 12 NJ tree" \
   -d "All pairs Distances" \
   -o /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/plots/chr12_LDCs_subsampled_only_ --no_labels \
   --highlight_samples /private/groups/patenlab/mira/centrolign/analysis/low_divergence_clades/max_dist_0.8_min_pairwise_0.95/chr12/subsampled_150.txt
```

![LDCs for chr12](../plots/chr12_LDCs_subsampled_only_pairwise_tree_heatmap.png)

Now, run centrolign on this subsampled set.

```sh
# Add 100kb to flanks of per sample alpha sat arrays

Regenerate fasta files with 100kb flanks

```sh
# Add 100 kb to bed files
cat /private/groups/patenlab/mira/centrolign/analysis/low_divergence_clades/max_dist_0.8_min_pairwise_0.95/chr12/subsampled_150.txt | cut -f2 -d"," | while read line ; do
  grep chr12 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/${line}_asat_arrays.bed | awk 'BEGIN{OFS="\t"} {$2=($2-100000<0)?0:$2-100000; $3=$3+100000; print}' > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds_100kb_flank/chr12/${line}_asat_arrays.100kb_flank.bed
done

# Prepare sample file
cat /private/groups/patenlab/mira/centrolign/analysis/low_divergence_clades/max_dist_0.8_min_pairwise_0.95/chr12/subsampled_150.txt | cut -f2 -d"," | sed 's/\./,/g' | while read line ; do grep $line /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_samples_with_asats.txt ; done >  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds_100kb_flank/chr12/chr12_subsample150.samples.txt

# Extract fasta files with flanks
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas_100kb_flank

mkdir -p logs
sbatch \
  /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/extract_fasta_r2_QCv2_100kb_flank_chr12_subsample_150.sh
```
3. Combine fastas for MSA

```sh
# check none are empty
find . -type f -empty

# combine 100kb fastas
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas_100kb_flank/chr12

ls | grep fasta | while read line ; do cat $line ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/runtime_debugging/r2_QC_v2_flanks_test/chr12/HPRC_r2_QC_v2_chr12_shuf50_100kb_flanks.fasta

# combine fastas
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas_100kb_flank/chr12

ls | while read line ; do
    cat $line
done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA_100kb_flank/chr12_subsample150/HPRC_r2_QC_v2_chr12_subsample150.fasta
```
4. Run centrolign MSA on same node

```sh
#!/bin/bash
#SBATCH --job-name=runtime_test_chr12_100kb_flanks
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=500gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodelist=phoenix-08
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=7-00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA_100kb_flank/chr12_subsample150/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr12_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA_100kb_flank/chr12_subsample150/HPRC_r2_QC_v2_chr12_subsample150.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA_100kb_flank/chr12_subsample150/chr12.subsample150.100kb_flanks.gfa
```

```
### Run Add Dummy Caps on all centrolign graphs for Glenn

```sh
#!/bin/bash
#SBATCH --job-name=add_dummy_caps
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=56gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=[1-24]%24
#SBATCH --exclude=phoenix-[09,10,22,23,24,18]
#SBATCH --output=logs/array_job_%A_task_%a.log
#SBATCH --time=1:00:00

CSV_FILE=/private/groups/patenlab/mira/centrolign/rGFA_tests/per_chrom_gfa.csv

CHR=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$CSV_FILE" | cut -f1 -d",")
GFA_LIST=$(awk -F"," "NR==$SLURM_ARRAY_TASK_ID" "$CSV_FILE" | cut -f2- -d",")
GFAS=`echo $GFA_LIST | tr -d '[]"' | tr ',' ' '`
OUTPUT_FILE=/private/groups/patenlab/mira/centrolign/rGFA_tests/centrolign_gfa_HPRC_r2_QC2_add_dummy_caps/${CHR}.centrolign.HPRC_r2_QCv2.gfa

echo $CHR
echo $GFAS

python3 /private/groups/patenlab/mira/centrolign/github/centromere-haplotype-sampling-pipeline/helper_scripts/add_dummy_caps.py \
    $GFAS -o $OUTPUT_FILE
```

### Extract Asat reads for all chromosomes

```sh
cd /private/groups/patenlab/mira/centrolign/rGFA_tests/extract_asat_reads
mkdir -p /private/groups/patenlab/mira/centrolign/rGFA_tests/extract_asat_reads/logs/

sbatch /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/extract_asat_reads.sh
```

Assessing file storage reqs:
```
grep -v -e 'chr4$' -e 'chr6$' -e 'chr9$' -e 'chr10$' -e 'chr11$' -e 'chr12$' -e 'chr17$' \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/*.bed \
    | wc -l

```
Resubmitting jobs that failed on phoenix 08
```
grep -l "Input/output error" /private/groups/patenlab/mira/centrolign/rGFA_tests/extract_asat_reads/logs/array_job_*.log \
    | sed 's/.*_task_\([0-9]*\)\.log/\1/' \
    | sort -n \
    | awk '
        NR==1 { start=$1; prev=$1; next }
        $1 == prev+1 { prev=$1; next }
        { if (start==prev) printf "%s,",start; else printf "%s-%s,",start,prev; start=$1; prev=$1 }
        END { if (start==prev) printf "%s\n",start; else printf "%s-%s\n",start,prev }
    '
7-11,13-19,21-58,60-76,78-155,217-232
```