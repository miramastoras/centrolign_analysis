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
