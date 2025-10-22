## Pairwise consistency benchmarking

The purpose of this analysis is to compare the consistency of the direct pairwise alignments and the pairwise alignments induced by the MSA runs.

### HPRC release 2 QC v2 - no flanks

#### Chr 5

```sh
#!/bin/bash
#SBATCH --job-name=chr5_subgroup0_pairwise_consistency
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign_%x.%j.log
#SBATCH --time=7-00:00

# subgroup 0
time python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/pairwise_consistency.py \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr5/subgroup_0/induced_pairwise_cigars/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr5/pairwise_cigar/pairwise_cigar_ \
    > /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/HPRC_r2_QCv2_chr5_subgroup0_pairwise_consistency.txt
```

```sh
#!/bin/bash
#SBATCH --job-name=chr5_subgroup1_pairwise_consistency
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign_%x.%j.log
#SBATCH --time=7-00:00

# subgroup 1
time python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/pairwise_consistency.py \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr5/subgroup_1/induced_pairwise_cigars/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr5/pairwise_cigar/pairwise_cigar_ \
    > /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/HPRC_r2_QCv2_chr5_subgroup1_pairwise_consistency.txt
```

#### Chr 12
```sh
#!/bin/bash
#SBATCH --job-name=chr12_subgroup1_pairwise_consistency
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign_%x.%j.log
#SBATCH --time=7-00:00

# subgroup 1
time python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/pairwise_consistency.py \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/induced_pairwise_cigars/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr12/pairwise_cigar/pairwise_cigar_ \
    > /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/HPRC_r2_QCv2_chr12_subgroup1_pairwise_consistency.txt
```
#### Chr 18

```sh
#!/bin/bash
#SBATCH --job-name=chr18_pairwise_consistency
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign_%x.%j.log
#SBATCH --time=7-00:00

# subgroup 1
time python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/pairwise_consistency.py \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr18/induced_pairwise_cigars/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr18/pairwise_cigar/pairwise_cigar_ \
    > /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/HPRC_r2_QCv2_chr18_pairwise_consistency.txt
```

#### Chr 17

```sh
#!/bin/bash
#SBATCH --job-name=chr17_pairwise_consistency
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign_%x.%j.log
#SBATCH --time=7-00:00

# subgroup 1
time python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/pairwise_consistency.py \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr17/induced_pairwise_cigars/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr17/pairwise_cigar/pairwise_cigar_ \
    > /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/HPRC_r2_QCv2_chr17_pairwise_consistency.txt
```

#### Chr 6
```sh
#!/bin/bash
#SBATCH --job-name=chr6_pairwise_consistency
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign_%x.%j.log
#SBATCH --time=7-00:00

# subgroup 1
time python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/pairwise_consistency.py \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr6/induced_pairwise_cigars/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr6/pairwise_cigar/pairwise_cigar_ \
    > /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/HPRC_r2_QCv2_chr6_pairwise_consistency.txt
```

#### Chr 12 subgroup 0
```sh
#!/bin/bash
#SBATCH --job-name=chr12_subgroup0_pairwise_consistency
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign_%x.%j.log
#SBATCH --time=7-00:00

# subgroup 0
time python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/pairwise_consistency.py \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_0/induced_pairwise_cigars/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr12/pairwise_cigar/pairwise_cigar_ \
    > /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/HPRC_r2_QCv2_chr12_subgroup0_pairwise_consistency.txt
```
#### Chr 8 subgroup 0

```sh
#!/bin/bash
#SBATCH --job-name=chr8_subgroup0_pairwise_consistency
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign_%x.%j.log
#SBATCH --time=7-00:00

# subgroup 0
time python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/pairwise_consistency.py \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr8/subgroup_0/induced_pairwise_cigars/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr8/pairwise_cigar/pairwise_cigar_ \
    > /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/HPRC_r2_QCv2_chr8_subgroup0_pairwise_consistency.txt
```
###

```sh
#!/bin/bash
#SBATCH --job-name=chr8_subgroup0_pairwise_consistency
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign_%x.%j.log
#SBATCH --time=7-00:00

# subgroup 0
time python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/pairwise_consistency.py \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr8/subgroup_0/induced_pairwise_cigars/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr8/pairwise_cigar/pairwise_cigar_ \
    > /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/HPRC_r2_QCv2_chr8_subgroup0_pairwise_consistency.txt
```
### Chr 4

```sh
#!/bin/bash
#SBATCH --job-name=chr4_pairwise_consistency
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign_%x.%j.log
#SBATCH --time=7-00:00


time python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/pairwise_consistency.py \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr4/induced_pairwise_cigars/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr4/pairwise_cigar/pairwise_cigar_ \
    > /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/HPRC_r2_QCv2_chr4_pairwise_consistency.txt
```
### Chr 7

Subgroup 0
```sh
#!/bin/bash
#SBATCH --job-name=chr7_subgroup0_pairwise_consistency
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign_%x.%j.log
#SBATCH --time=7-00:00


time python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/pairwise_consistency.py \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr7/subgroup_0/induced_pairwise_cigars/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr7/pairwise_cigar/pairwise_cigar_ \
    > /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/HPRC_r2_QCv2_chr7_subgroup0_pairwise_consistency.txt

```
Subgroup 1

```sh
#!/bin/bash
#SBATCH --job-name=chr7_subgroup1_pairwise_consistency
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign_%x.%j.log
#SBATCH --time=7-00:00


time python3 /private/groups/patenlab/mira/centrolign/github/centromere-scripts/benchmarking/pairwise_consistency.py \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr7/subgroup_1/induced_pairwise_cigars/pairwise_cigar_ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr7/pairwise_cigar/pairwise_cigar_ \
    > /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/HPRC_r2_QCv2_chr7_subgroup1_pairwise_consistency.txt
```

### Get patristic distances from the input NJ tree

```sh
chromosomes=("chr4" "chr5" "chr6" "chr8" "chr7" "chr12" "chr17" "chr18")

for chr in "${chromosomes[@]}"; do
  /private/home/mmastora/progs/centrolign/build/tree_pair_dist /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/${chr}_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk > /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/tree_pair_dist/${chr}_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk.pair_dists.tsv
done
```

### Pairwise consistency per-chrom histograms

```sh
chromosomes=("chr12" "chr5" "chr6")

for chr in "${chromosomes[@]}"; do
  Rscript /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/self_consistency.R \
    /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/HPRC_r2_QCv2_${chr}_pairwise_consistency.txt \
    /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/tree_pair_dist/${chr}_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk.pair_dists.tsv \
    ${chr} \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/${chr}_r2_QC_v2_centrolign_pairwise_distance.csv \
    /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/plots/HPRC_r2_QCv2
done
```

#### Synteny plots for sanity check

cd /private/groups/patenlab/mira/centrolign/analysis/chr3_chr4_interspersed_hsat_r2

synteny plot
```sh
cd /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/chr6_examples

samples=("HG00290.1" "HG06807.1")

# Process each line in the fai file
for smp in "${samples[@]}"; do
  samtools faidx /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/chr6/${smp}_chr6_hor_array.fasta
  while read -r seqid length rest; do
      # Create a BED file for each sequence
      # BED format: chrom start end name score strand thickStart thickEnd itemRgb
      echo -e "${seqid}\t0\t${length}\tentry\t0\t.\t0\t${length}\t0,0,0" > dummy_beds/"${seqid}.bed"
      echo "Created bed file for ${seqid}"
  done < /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/chr6/${smp}_chr6_hor_array.fasta.fai
done


python /private/groups/migalab/juklucas/centrolign/chr12_test125/synteny_plot_bokeh.py   \
    --beds \
      /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/chr6_examples/dummy_beds/HG00290.1.bed \
      /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/chr6_examples/dummy_beds/HG06807.1.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/chr6/pairwise_cigar/pairwise_cigar_HG00290.1_HG06807.1.txt \
    --output /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/chr6_examples/plots/synteny_direct_HG00290.1_HG06807.1.html \
    --web

#!/bin/sh
python /private/groups/migalab/juklucas/centrolign/chr12_test125/synteny_plot_bokeh.py   \
    --beds \
      /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/chr6_examples/dummy_beds/HG00290.1.bed \
      /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/chr6_examples/dummy_beds/HG06807.1.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr6/induced_pairwise_cigars/pairwise_cigar_HG00290.1_HG06807.1.txt \
    --output /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/chr6_examples/plots/synteny_induced_HG00290.1_HG06807.1.html \
    --web
```
