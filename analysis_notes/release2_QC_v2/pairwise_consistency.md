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
    > /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/HPRC_r2_QCv2_chr18_subgroup1_pairwise_consistency.txt
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
    > /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/HPRC_r2_QCv2_chr6_subgroup1_pairwise_consistency.txt
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
