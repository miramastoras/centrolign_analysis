## Submitting centrolign MSA for HPRC release 2 QC v2, using all pairs neighbor-joining trees split manually for runtime

### Chr 5

Description of subgroups: https://docs.google.com/presentation/d/1V6O15S6s3WA9YeKI7O4qbRdZJjpjFO4xSXR5j8g1dCs/edit?slide=id.g38e95b68347_0_48#slide=id.g38e95b68347_0_48


#### Subgroup 0:

Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr5/subgroup_0/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr5/subgroup_0/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr5/subgroup_0/

```
Slurm script:

```sh
#!/bin/bash
#SBATCH --job-name=chr5_MSA_subgroup_0
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=800gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=7-00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr5/subgroup_0/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr5/split_by_2/subgroup_0_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr5/split_by_2/subgroup_0_seqs.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr5/subgroup_0/chr5.subgroup_0.centrolign.gfa
```

#### Subgroup 1:

Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr5/subgroup_1/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr5/subgroup_1/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr5/subgroup_1/

```
Slurm script:

```sh
#!/bin/bash
#SBATCH --job-name=chr5_MSA_subgroup_1
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=800gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=7-00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr5/subgroup_1/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr5/split_by_2/subgroup_1_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr5/split_by_2/subgroup_1_seqs.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr5/subgroup_1/chr5.subgroup_1.centrolign.gfa
```
### Chr 8

Description of subgroups: https://docs.google.com/presentation/d/1V6O15S6s3WA9YeKI7O4qbRdZJjpjFO4xSXR5j8g1dCs/edit?slide=id.g38e95b68347_0_48#slide=id.g38e95b68347_0_48

#### Subgroup 0:

Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr8/subgroup_0/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr8/subgroup_0/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr8/subgroup_0/

```
Slurm script:

```sh
#!/bin/bash
#SBATCH --job-name=chr8_MSA_subgroup_0
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=800gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=7-00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr8/subgroup_0/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr8/split_by_2/subgroup_0_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr8/split_by_2/subgroup_0_seqs.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr8/subgroup_0/chr8.subgroup_0.centrolign.gfa
```

#### Subgroup 1:

Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr8/subgroup_1/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr8/subgroup_1/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr8/subgroup_1/

```
Slurm script:

```sh
#!/bin/bash
#SBATCH --job-name=chr8_MSA_subgroup_1
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=800gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=7-00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr8/subgroup_1/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr8/split_by_2/subgroup_1_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr8/split_by_2/subgroup_1_seqs.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr8/subgroup_1/chr8.subgroup_1.centrolign.gfa
```
## Chr 12

#### Subgroup 0:

Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_0/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_0/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_0/

```
Slurm script:

```sh
#!/bin/bash
#SBATCH --job-name=chr12_MSA_subgroup_0
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=800gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=7-00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_0/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr12/split_by_2/subgroup_0_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr12/split_by_2/subgroup_0_seqs.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_0/chr12.subgroup_0.centrolign.gfa
```
#### Subgroup 1:

Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/

```
Slurm script:

```sh
#!/bin/bash
#SBATCH --job-name=chr12_MSA_subgroup_1
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=800gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=7-00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr12/split_by_2/subgroup_1_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr12/split_by_2/subgroup_1_seqs.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/chr12.subgroup_1.centrolign.gfa
```
