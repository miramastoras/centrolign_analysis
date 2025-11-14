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
Inducing pairwise alignments
```sh
#!/bin/bash
#SBATCH --job-name=chr5_MSA_subgroup_0_induced_pairwise
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=400gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=7-00:00

time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr5/subgroup_0/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr5/split_by_2/subgroup_0_tree.nwk \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr5/subgroup_0/induced_pairwise_cigars/pairwise_cigar \
  -R \
  --threads 32 \
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
#SBATCH --job-name=chr5_MSA_subgroup_1_induced_pairwise
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=800gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=7-00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr5/subgroup_1/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr5/split_by_2/subgroup_1_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr5/split_by_2/subgroup_1_seqs.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr5/subgroup_1/chr5.subgroup_1.centrolign.gfa
```

Induced pairwise cigars
```sh
#!/bin/bash
#SBATCH --job-name=chr5_MSA_subgroup_1_induced_pairwise
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=400gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=7-00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr5/subgroup_1/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr5/split_by_2/subgroup_1_tree.nwk \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr5/subgroup_1/induced_pairwise_cigars/pairwise_cigar \
  -R \
  --threads 32 \
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
Induced pairwise cigars

```sh
#!/bin/bash
#SBATCH --job-name=chr8_MSA_subgroup_0_induced_pairwise
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/induced_pairwise_centrolign_%x.%j.log
#SBATCH --time=12:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr8/subgroup_0/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr8/split_by_2/subgroup_0_tree.nwk \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr8/subgroup_0/induced_pairwise_cigars/pairwise_cigar \
  -R \
  --threads 32 \
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
Induced pairwise cigars

```sh
#!/bin/bash
#SBATCH --job-name=chr8_MSA_subgroup_1_induced_pairwise

#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=12:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr8/subgroup_1/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr8/split_by_2/subgroup_1_tree.nwk \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr8/subgroup_1/induced_pairwise_cigars/pairwise_cigar \
  -R \
  --threads 32 \
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
Pairwise cigars
```sh
#!/bin/bash
#SBATCH --job-name=chr12_MSA_subgroup_0_induced_pairwise
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=12:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_0/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr12/split_by_2/subgroup_0_tree.nwk \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_0/induced_pairwise_cigars/pairwise_cigar \
  -R \
  --threads 32 \
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
Induced pairwise cigars

```sh
#!/bin/bash
#SBATCH --job-name=chr12_MSA_subgroup_1_induced_pairwise
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=400gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=7-00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr12/split_by_2/subgroup_1_tree.nwk \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/induced_pairwise_cigars/pairwise_cigar \
  -R \
  --threads 32 \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr12/split_by_2/subgroup_1_seqs.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/chr12.subgroup_1.centrolign.gfa
```

### Chr17, Chr 18, Chr 6 - small sample size so no splitting needed

Combine fasta file for centrolign
```sh
chromosomes=("chr17" "chr18" "chr6")

for chr in "${chromosomes[@]}"
do
  cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/fasta_lists/release2_QC_v2_${chr}.txt | while read line ; do cat $line ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.${chr}.fasta

  samtools faidx /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.${chr}.fasta

done
```
##### Chr17

Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr17/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr17/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr17/
```

```sh
#!/bin/bash
#SBATCH --job-name=chr17_MSA
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr17/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr17_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.chr17.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr17/chr17.centrolign.gfa
```
Induced pairwise alignments

```sh
#!/bin/bash
#SBATCH --job-name=chr17_MSA_induced_pairwise
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/induced_pairwise_centrolign_%x.%j.log
#SBATCH --time=12:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr17/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr17_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr17/induced_pairwise_cigars/pairwise_cigar \
  -R \
  --threads 32 \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.chr17.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr17/chr17.centrolign.gfa
```

##### Chr18

Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr18/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr18/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr18/
```

```sh
#!/bin/bash
#SBATCH --job-name=chr18_MSA
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr18/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr18_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.chr18.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr18/chr18.centrolign.gfa
```

Induced pairwise cigars
```sh
#!/bin/bash
#SBATCH --job-name=chr18_MSA_induced_pairwise
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=300gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/induced_pairwise_centrolign_%x.%j.log
#SBATCH --time=12:00:00

time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr18/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr18_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr18/induced_pairwise_cigars/pairwise_cigar \
  -R \
  --threads 32 \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.chr18.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr18/chr18.centrolign.gfa
```
##### Chr6

Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr6/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr6/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr6/
```

```sh
#!/bin/bash
#SBATCH --job-name=chr6_MSA
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr6/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr6_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.chr6.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr6/chr6.centrolign.gfa
```
Induced pairwise cigars
```sh
#!/bin/bash
#SBATCH --job-name=chr6_MSA_induced_pairwise
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/induced_pairwise_centrolign_%x.%j.log
#SBATCH --time=12:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr6/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr6_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.chr6.fasta \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr6/induced_pairwise_cigars/pairwise_cigar \
  -R \
  --threads 32 \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr6/chr6.centrolign.gfa
```


##### Chr2 - Submitting Chr 2 as a single GFA, because we know we want to keep these together

Combine all fastas
```sh
chr=chr2

cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/fasta_lists/release2_QC_v2_${chr}.txt | while read line ; do cat $line ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.${chr}.fasta

samtools faidx /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.${chr}.fasta
```
Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr2/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr2/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr2/
```

```sh
#!/bin/bash
#SBATCH --job-name=chr2_MSA
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr2/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr2_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.chr2.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr2/chr2.centrolign.gfa
```

#### Chr 4

Combine all fastas
```sh
chr=chr4

cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/fasta_lists/release2_QC_v2_${chr}.txt | while read line ; do cat $line ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.${chr}.fasta

samtools faidx /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.${chr}.fasta
```
Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr4/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr4/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr4/
```

```sh
#!/bin/bash
#SBATCH --job-name=chr4_MSA
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr4/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr4_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.chr4.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr4/chr4.centrolign.gfa
```
Induced pairwise cigars
```sh
#!/bin/bash
#SBATCH --job-name=chr4_MSA_induced_pairwise
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/centrolign_induced_pairwise_%x.%j.log
#SBATCH --time=12:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr4/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr4_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr4/induced_pairwise_cigars/pairwise_cigar \
  -R \
  --threads 32 \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.chr4.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr4/chr4.centrolign.gfa
```

#### Chr 3

Combine all fastas
```sh
chr=chr3

cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/fasta_lists/release2_QC_v2_${chr}.txt | while read line ; do cat $line ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.${chr}.fasta

samtools faidx /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.${chr}.fasta
```
Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr3/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr3/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr3/
```

```sh
#!/bin/bash
#SBATCH --job-name=chr3_MSA
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr3/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr3_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.chr3.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr3/chr3.centrolign.gfa
```

Induced pairwise cigars
```sh
#!/bin/bash
#SBATCH --job-name=chr3_MSA_induced_pairwise
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=12:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr3/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr3_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr3/induced_pairwise_cigars/pairwise_cigar \
  -R \
  --threads 32 \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.chr3.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr3/chr3.centrolign.gfa
```

#### Chr X - running all together

Combine all fastas
```sh
chr=chrX

cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/fasta_lists/release2_QC_v2_${chr}.txt | while read line ; do cat $line ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.${chr}.fasta

samtools faidx /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.${chr}.fasta
```
Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrX/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrX/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrX/
```

```sh
#!/bin/bash
#SBATCH --job-name=chrX_MSA
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrX/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chrX_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.chrX.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrX/chrX.centrolign.gfa
```
#### Chr 7

#### Subgroup 0:

Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr7/subgroup_0/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr7/subgroup_0/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr7/subgroup_0/

```
Slurm script:

```sh
#!/bin/bash
#SBATCH --job-name=chr7_MSA_subgroup_0
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr7/subgroup_0/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr7/split_by_2/subgroup_0_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr7/split_by_2/subgroup_0_seqs.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr7/subgroup_0/chr7.subgroup_0.centrolign.gfa
```
Induced pairwise alignments
```sh
#!/bin/bash
#SBATCH --job-name=chr7_MSA_subgroup_0_induced_pairwise
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=12:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr7/subgroup_0/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr7/split_by_2/subgroup_0_tree.nwk \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr7/subgroup_0/induced_pairwise_cigars/pairwise_cigar \
  -R \
  --threads 32 \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr7/split_by_2/subgroup_0_seqs.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr7/subgroup_0/chr7.subgroup_0.centrolign.gfa
```

#### Subgroup 1:

Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr7/subgroup_1/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr7/subgroup_1/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr7/subgroup_1/

```
Slurm script:

```sh
#!/bin/bash
#SBATCH --job-name=chr7_MSA_subgroup_1
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr7/subgroup_1/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr7/split_by_2/subgroup_1_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr7/split_by_2/subgroup_1_seqs.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr7/subgroup_1/chr7.subgroup_1.centrolign.gfa
```
Induced pairwise

```sh
#!/bin/bash
#SBATCH --job-name=chr7_MSA_subgroup_1_induced_pairwise
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=7-00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr7/subgroup_1/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr7/split_by_2/subgroup_1_tree.nwk \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr7/subgroup_1/induced_pairwise_cigars/pairwise_cigar \
  -R \
  --threads 32 \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr7/split_by_2/subgroup_1_seqs.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr7/subgroup_1/chr7.subgroup_1.centrolign.gfa
```

#### Chr 22

#### Subgroup 0:

Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr22/subgroup_0/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr22/subgroup_0/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr22/subgroup_0/

```
Slurm script:

```sh
#!/bin/bash
#SBATCH --job-name=chr22_MSA_subgroup_0
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr22/subgroup_0/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr22/split_by_2/subgroup_0_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr22/split_by_2/subgroup_0_seqs.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr22/subgroup_0/chr22.subgroup_0.centrolign.gfa
```
Induced pairwise

```sh
#!/bin/bash
#SBATCH --job-name=chr22_MSA_subgroup_0_induced_pairwise
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=12:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr22/subgroup_0/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr22/split_by_2/subgroup_0_tree.nwk \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr22/subgroup_0/induced_pairwise_cigars/pairwise_cigar \
  -R \
  --threads 32 \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr22/split_by_2/subgroup_0_seqs.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr22/subgroup_0/chr22.subgroup_0.centrolign.gfa
```


#### Subgroup 1:

Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr22/subgroup_1/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr22/subgroup_1/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr22/subgroup_1/

```
Slurm script:

```sh
#!/bin/bash
#SBATCH --job-name=chr22_MSA_subgroup_1
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr22/subgroup_1/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr22/split_by_2/subgroup_1_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr22/split_by_2/subgroup_1_seqs.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr22/subgroup_1/chr22.subgroup_1.centrolign.gfa
```

induced pairwise cigars

```sh
#!/bin/bash
#SBATCH --job-name=chr22_MSA_subgroup_1_induced_pairwise
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=12:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr22/subgroup_1/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr22/split_by_2/subgroup_1_tree.nwk \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr22/subgroup_1/induced_pairwise_cigars/pairwise_cigar \
  -R \
  --threads 32 \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr22/split_by_2/subgroup_1_seqs.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr22/subgroup_1/chr22.subgroup_1.centrolign.gfa
```



### Chr 20

Decided to run entire sample list together, because of how fast chr7 ran

Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr20/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr20/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr20/
```

```sh
#!/bin/bash
#SBATCH --job-name=chr20_MSA
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr20/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr20_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.chr20.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr20/chr20.centrolign.gfa
```
Running the whole thing together was a bad idea - out of memory after 1500 Gb!

Checking the last two subproblems to see if they contain all samples
```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr20/subproblems
tail -n 1 _info.txt | cut -f2 | sed 's/,/\n/g' | wc -l #180
tail -n 1 _info.txt | cut -f2 | sed 's/,/\n/g' > /private/groups/patenlab/mira/_D7FF81E954CBC598.gfa.samples.txt

tail -n 2 _info.txt  | head -n 1 | cut -f2 | sed 's/,/\n/g' | wc -l #149
tail -n 2 _info.txt  | head -n 1 | cut -f2 | sed 's/,/\n/g' > /private/groups/patenlab/mira/_B68D67B1DDE9835E.gfa.samples.txt
```
Plotting the heatmaps for these samples [heatmaps](https://docs.google.com/presentation/d/1V6O15S6s3WA9YeKI7O4qbRdZJjpjFO4xSXR5j8g1dCs/edit?slide=id.g3a0bfb799c2_0_0#slide=id.g3a0bfb799c2_0_0)

```sh
chromosomes=("chr20")

for chr in "${chromosomes[@]}"
do
  python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap_v2.py \
    -t /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_all_pairs_nj_tree.format5.nwk \
    -s /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}.samples.txt \
    -p /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_pairwise_distance.csv \
    -m "Centrolign all pairs distances" \
    -n "${chr} NJ tree" \
    -d "All pairs Distances" \
    -o /Users/miramastoras/Desktop/color_subgroups_heatmap/${chr}__D7FF81E954CBC598gfa --no_labels \
    --highlight_samples /Users/miramastoras/Desktop/color_subgroups_heatmap/_D7FF81E954CBC598.gfa.samples.txt
done

chromosomes=("chr20")

for chr in "${chromosomes[@]}"
do
  python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap_v2.py \
    -t /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_all_pairs_nj_tree.format5.nwk \
    -s /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}.samples.txt \
    -p /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_pairwise_distance.csv \
    -m "Centrolign all pairs distances" \
    -n "${chr} NJ tree" \
    -d "All pairs Distances" \
    -o /Users/miramastoras/Desktop/color_subgroups_heatmap/${chr}_B68D67B1DDE9835Egfa --no_labels \
    --highlight_samples /Users/miramastoras/Desktop/color_subgroups_heatmap/_B68D67B1DDE9835E.gfa.samples.txt
done
```

Confirmed that all samples are contained in these two subproblems!

Renaming for clarity:
```sh
# SUBGROUP D - 180 samples
/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr20/subproblems/_D7FF81E954CBC598.gfa

# SUBGROUP B - 149 samples
/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr20/subproblems/_B68D67B1DDE9835E.gfa
```

Getting the induced pairwise cigar strings from these subproblems:
```sh
# get sample lists
tail -n 1 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr20/subproblems/_info.txt | cut -f2 | sed 's/,/\n/g' > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr20/induced_pairwise_cigars/chr20_subgroup_D.samples.txt

tail -n 2 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr20/subproblems/_info.txt | head -n 1 | cut -f2 | sed 's/,/\n/g' > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr20/induced_pairwise_cigars/chr20_subgroup_B.samples.txt

# create new fasta with just those samples

## subgroup D
cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr20/induced_pairwise_cigars/chr20_subgroup_D.samples.txt | while read line ; do
    cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/chr20/${line}_chr20_hor_array.fasta >> /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr20/induced_pairwise_cigars/chr20_subgroup_D.fasta
done

## subgroup B
cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr20/induced_pairwise_cigars/chr20_subgroup_B.samples.txt | while read line ; do
    cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/chr20/${line}_chr20_hor_array.fasta >> /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr20/induced_pairwise_cigars/chr20_subgroup_B.fasta
done
```

```sh
#!/bin/bash
#SBATCH --job-name=chr20_MSA_subgroupB
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=12:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr20/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr20_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr20/induced_pairwise_cigars/subgroup_B/pairwise_cigar \
  -R \
  --threads 32 \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr20/induced_pairwise_cigars/chr20_subgroup_B.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr20/chr20.centrolign.gfa
```

```sh
#!/bin/bash
#SBATCH --job-name=chr20_MSA_subgroupD
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=12:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr20/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr20_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr20/induced_pairwise_cigars/subgroup_D/pairwise_cigar \
  -R \
  --threads 32 \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr20/induced_pairwise_cigars/chr20_subgroup_D.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr20/chr20.centrolign.gfa
```

### Chr 1

Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr1/subgroup_A/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr1/subgroup_A/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr1/subgroup_A
```

```sh
#!/bin/bash
#SBATCH --job-name=chr1_subgroup_A_MSA
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr1/subgroup_A/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr1_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr1/combine_final_subgroups/chr1.subgroup_A.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr1/subgroup_A/chr1_subgroup_A.centrolign.gfa
```
Chr 1 subgroup A maxed out at 1400 Gb of memory.

Checking the last subproblem to see which samples it contains
```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr1/subgroup_A/subproblems
tail -n 1 _info.txt | cut -f2 | sed 's/,/\n/g' | wc -l #100
tail -n 1 _info.txt | cut -f2 | sed 's/,/\n/g' > /private/groups/patenlab/mira/color_subgroups_heatmap/_598707B3F6931A5A.gfa.samples.txt
```
Plotting on heatmap
```sh
chromosomes=chr1
for chr in "${chromosomes[@]}"
do
  python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap_v2.py \
    -t /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_all_pairs_nj_tree.format5.nwk \
    -s /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}.samples.txt \
    -p /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_pairwise_distance.csv \
    -m "Centrolign all pairs distances" \
    -n "${chr} NJ tree" \
    -d "All pairs Distances" \
    -o /Users/miramastoras/Desktop/github_repos/centrolign_analysis/analysis_notes/release2_QC_v2/plots/${chr}_r2_QC_v2_subgroupA_598707B3F6931A5A.gfa_ --no_labels \
    --highlight_samples /Users/miramastoras/Desktop/color_subgroups_heatmap/_598707B3F6931A5A.gfa.samples.txt
done
```

Now checking the next subgroup in the info file that has seemingly non overlapping samples
```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr1/subgroup_A/subproblems
tail -n 1 _info.txt | cut -f2 | sed 's/,/\n/g' | wc -l #100
tail -n 1 _info.txt | cut -f2 | sed 's/,/\n/g' > /private/groups/patenlab/mira/color_subgroups_heatmap/_598707B3F6931A5A.gfa.samples.txt

grep /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr1/subgroup_A/subproblems/_476C7B62F401E99B.gfa _info.txt | cut -f2 | sed 's/,/\n/g' | wc -l #79

grep /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr1/subgroup_A/subproblems/_476C7B62F401E99B.gfa _info.txt | cut -f2 | sed 's/,/\n/g' | sed 's/,/\n/g' > /private/groups/patenlab/mira/color_subgroups_heatmap/_476C7B62F401E99B.gfa.samples.txt
```
Plotting on heatmap
```sh
chromosomes=chr1
for chr in "${chromosomes[@]}"
do
  python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap_v2.py \
    -t /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_all_pairs_nj_tree.format5.nwk \
    -s /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}.samples.txt \
    -p /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_pairwise_distance.csv \
    -m "Centrolign all pairs distances" \
    -n "${chr} NJ tree" \
    -d "All pairs Distances" \
    -o /Users/miramastoras/Desktop/github_repos/centrolign_analysis/analysis_notes/release2_QC_v2/plots/${chr}_r2_QC_v2_subgroupA_476C7B62F401E99B.gfa_ --no_labels \
    --highlight_samples /Users/miramastoras/Desktop/color_subgroups_heatmap/_476C7B62F401E99B.gfa.samples.txt
done
```


Subgroup B

Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr1/subgroup_B/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr1/subgroup_B/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr1/subgroup_B
```

```sh
#!/bin/bash
#SBATCH --job-name=chr1_subgroup_B_MSA
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr1/subgroup_B/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr1_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr1/combine_final_subgroups/chr1.subgroup_B.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr1/subgroup_B/chr1_subgroup_B.centrolign.gfa
```

Induced pairwise

```sh
#!/bin/bash
#SBATCH --job-name=chr1_subgroup_B_MSA_induced_pairwise
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=12:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr1/subgroup_B/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr1_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr1/subgroup_B/induced_pairwise_cigars/pairwise_cigar \
  -R \
  --threads 32 \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr1/combine_final_subgroups/chr1.subgroup_B.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr1/subgroup_B/chr1_subgroup_B.centrolign.gfa
```

### Chr Y

Combine all fastas
```sh
chr=chrY

cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/fasta_lists/release2_QC_v2_${chr}.txt | while read line ; do cat $line ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.${chr}.fasta

samtools faidx /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.${chr}.fasta
```
Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrY/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrY/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrY/
```


```sh
#!/bin/bash
#SBATCH --job-name=chrY_MSA
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrY/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chrY_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.chrY.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrY/chrY.centrolign.gfa
```

ChrY
```sh
#!/bin/bash
#SBATCH --job-name=chrY_MSA_induced_pairwise
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=12:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrY/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chrY_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrY/induced_pairwise_cigars/pairwise_cigar \
  -R \
  --threads 32 \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/full_fastas/HPRC_release2_QC_v2.chrY.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrY/chrY.centrolign.gfa
```

### Chr 19

Subgroup 0

Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr19/subgroup_0/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr19/subgroup_0/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr19/subgroup_0/

```
Slurm script:

```sh
#!/bin/bash
#SBATCH --job-name=chr19_MSA_subgroup_0
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr19/subgroup_0/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr19/split_by_2/subgroup_0_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr19/split_by_2/subgroup_0_seqs.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr19/subgroup_0/chr19.subgroup_0.centrolign.gfa
```

Induced pairwise cigars
```sh
#!/bin/bash
#SBATCH --job-name=chr19_MSA_subgroup_0_induced_pairwise
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=12:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr19/subgroup_0/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr19/split_by_2/subgroup_0_tree.nwk \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr19/subgroup_0/induced_pairwise_cigars/pairwise_cigar \
  -R \
  --threads 32 \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr19/split_by_2/subgroup_0_seqs.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr19/subgroup_0/chr19.subgroup_0.centrolign.gfa
```

Subgroup 1

Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr19/subgroup_1/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr19/subgroup_1/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr19/subgroup_1/

```
Slurm script:

```sh
#!/bin/bash
#SBATCH --job-name=chr19_MSA_subgroup_1
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr19/subgroup_1/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr19/split_by_2/subgroup_1_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr19/split_by_2/subgroup_1_seqs.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr19/subgroup_1/chr19.subgroup_1.centrolign.gfa
```

Induced pairwise cigars

```sh
#!/bin/bash
#SBATCH --job-name=chr19_MSA_subgroup_1_induced_pairwise
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=12:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr19/subgroup_1/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr19/split_by_2/subgroup_1_tree.nwk \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr19/subgroup_1/induced_pairwise_cigars/pairwise_cigar \
  -R \
  --threads 32 \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr19/split_by_2/subgroup_1_seqs.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr19/subgroup_1/chr19.subgroup_1.centrolign.gfa
```

### Chr 21

Subgroup 0

Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr21/subgroup_0/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr21/subgroup_0/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr21/subgroup_0/

```
Slurm script:

```sh
#!/bin/bash
#SBATCH --job-name=chr21_MSA_subgroup_0
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr21/subgroup_0/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr21/split_by_2/subgroup_0_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr21/split_by_2/subgroup_0_seqs.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr21/subgroup_0/chr21.subgroup_0.centrolign.gfa
```
Induced pairwise

```sh
#!/bin/bash
#SBATCH --job-name=chr21_MSA_subgroup_0_induced_pairwise
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=12:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr21/subgroup_0/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr21/split_by_2/subgroup_0_tree.nwk \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr21/subgroup_0/induced_pairwise_cigars/pairwise_cigar \
  -R \
  --threads 32 \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr21/split_by_2/subgroup_0_seqs.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr21/subgroup_0/chr21.subgroup_0.centrolign.gfa
```

Subgroup 1

Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr21/subgroup_1/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr21/subgroup_1/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr21/subgroup_1/

```
Slurm script:

```sh
#!/bin/bash
#SBATCH --job-name=chr21_MSA_subgroup_1
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr21/subgroup_1/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr21/split_by_2/subgroup_1_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr21/split_by_2/subgroup_1_seqs.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr21/subgroup_1/chr21.subgroup_1.centrolign.gfa
```

```sh
#!/bin/bash
#SBATCH --job-name=chr21_MSA_subgroup_1_induced_pairwise
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=12:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr21/subgroup_1/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr21/split_by_2/subgroup_1_tree.nwk \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr21/subgroup_1/induced_pairwise_cigars/pairwise_cigar \
  -R \
  --threads 32 \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr21/split_by_2/subgroup_1_seqs.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr21/subgroup_1/chr21.subgroup_1.centrolign.gfa
```

Induced pairwise


### Chr10

Subgroup A
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr10/subgroup_A/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr10/subgroup_A/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr10/subgroup_A/

```
Slurm script
```sh
#!/bin/bash
#SBATCH --job-name=chr10_MSA_subgroup_A
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr10/subgroup_A/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr10_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr10/combine_final_subgroups/chr10.subgroup_A.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr10/subgroup_A/chr10.subgroup_A.centrolign.gfa
```

Induced pairwise cigars
```sh
#!/bin/bash
#SBATCH --job-name=chr10_MSA_subgroup_A_induced_pairwise
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=12:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr10/subgroup_A/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr10_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr10/subgroup_A/induced_pairwise_cigars/pairwise_cigar \
  -R \
  --threads 32 \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr10/combine_final_subgroups/chr10.subgroup_A.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr10/subgroup_A/chr10.subgroup_A.centrolign.gfa
```

Subgroup B
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr10/subgroup_B/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr10/subgroup_B/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr10/subgroup_B/

```
Slurm script
```sh
#!/bin/bash
#SBATCH --job-name=chr10_MSA_subgroup_B
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr10/subgroup_B/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr10_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr10/combine_final_subgroups/chr10.subgroup_B.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr10/subgroup_B/chr10.subgroup_B.centrolign.gfa
```
Induced Pairwise cigars
```sh
#!/bin/bash
#SBATCH --job-name=chr10_MSA_subgroup_B
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=12:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr10/subgroup_B/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr10_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr10/subgroup_B/induced_pairwise_cigars/pairwise_cigar \
  -R \
  --threads 32 \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr10/combine_final_subgroups/chr10.subgroup_B.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr10/subgroup_B/chr10.subgroup_B.centrolign.gfa

```

### Chr9


Subgroup 0

Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr9/subgroup_0/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr9/subgroup_0/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr9/subgroup_0/

```
Slurm script:

```sh
#!/bin/bash
#SBATCH --job-name=chr9_MSA_subgroup_0
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr9/subgroup_0/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr9/split_by_2/subgroup_0_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr9/split_by_2/subgroup_0_seqs.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr9/subgroup_0/chr9.subgroup_0.centrolign.gfa
```


Subgroup 1

Working directory:
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr9/subgroup_1/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr9/subgroup_1/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr9/subgroup_1/

```
Slurm script:

```sh
#!/bin/bash
#SBATCH --job-name=chr9_MSA_subgroup_1
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr9/subgroup_1/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr9/split_by_2/subgroup_1_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr9/split_by_2/subgroup_1_seqs.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr9/subgroup_1/chr9.subgroup_1.centrolign.gfa
```

induced pairwise cigars
```sh
#!/bin/bash
#SBATCH --job-name=ch9_MSA_subgroup_1_induced_pairwise
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=12:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr9/subgroup_1/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr9/split_by_2/subgroup_1_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr9/split_by_2/subgroup_1_seqs.fasta \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr9/subgroup_1/induced_pairwise_cigars/pairwise_cigar \
  -R \
  --threads 32 \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr9/subgroup_1/chr9.subgroup_1.centrolign.gfa
```

### Chr 11

Subgroup A
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr11/subgroup_A/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr11/subgroup_A/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr11/subgroup_A/

```
Slurm script
```sh
#!/bin/bash
#SBATCH --job-name=chr11_MSA_subgroup_A
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr11/subgroup_A/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr11_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr11/combine_final_subgroups/chr11.subgroup_A.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr11/subgroup_A/chr11.subgroup_A.centrolign.gfa
```

Subgroup B
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr11/subgroup_B/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr11/subgroup_B/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr11/subgroup_B/

```
Slurm script
```sh
#!/bin/bash
#SBATCH --job-name=chr11_MSA_subgroup_B
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr11/subgroup_B/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr11_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr11/combine_final_subgroups/chr11.subgroup_B.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr11/subgroup_B/chr11.subgroup_B.centrolign.gfa
```

Induced pairwise cigars
Slurm script
```sh
#!/bin/bash
#SBATCH --job-name=chr11_MSA_subgroup_B_induced_pairwise
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=12:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr11/subgroup_B/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr11_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr11/subgroup_B/induced_pairwise_cigars/pairwise_cigar \
  -R \
  --threads 32 \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr11/combine_final_subgroups/chr11.subgroup_B.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr11/subgroup_B/chr11.subgroup_B.centrolign.gfa
```


Subgroup C
```sh
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr11/subgroup_C/subproblems/
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr11/subgroup_C/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr11/subgroup_C/
```

Slurm script
```sh
#!/bin/bash
#SBATCH --job-name=chr11_MSA_subgroup_C
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr11/subgroup_C/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr11_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr11/combine_final_subgroups/chr11.subgroup_C.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr11/subgroup_C/chr11.subgroup_C.centrolign.gfa
```

Induced pairwise cigars

Slurm script
```sh
#!/bin/bash
#SBATCH --job-name=chr11_MSA_subgroup_C_induced_pairwise
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/centrolign_%x.%j.log
#SBATCH --time=12:00:00


time /private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr11/subgroup_C/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr11_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr11/subgroup_C/induced_pairwise_cigars/pairwise_cigar \
  -R \
  --threads 32 \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/split_nj_trees/chr11/combine_final_subgroups/chr11.subgroup_C.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr11/subgroup_C/chr11.subgroup_C.centrolign.gfa
```
### Backing up the MSAs

Copy to new folder then tar them up
```sh
mkdir -p /private/groups/patenlab/mira/centrolign_MSAs_r2_QC_v2/

tail -n 34 /private/groups/patenlab/mira/centrolign/analysis/pairwise_consistency/centrolign_results.csv | cut -f 9 -d"," | sort | tail -n 31 | while read line ; do
    cp $line /private/groups/patenlab/mira/centrolign_MSAs_r2_QC_v2/
  done

cd

tar -zcvf centrolign_MSAs_r2_QC_v2.tar.gz centrolign_MSAs_r2_QC_v2/
```

Uploaded to my personal google drive: https://drive.google.com/file/d/1YGObICA_RHaKWkR7UID7rY8GXZzC2FHf/view?usp=drive_link 
