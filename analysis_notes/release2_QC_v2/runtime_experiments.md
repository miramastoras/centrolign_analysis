### Experiments to test runtime of centrolign

#### Does adding 100kb of flank really increase the runtime drastically?

For chr12, select 50 samples. Run with 100 kb of flanks, and without 100 kb of flanks

1. Randomly select 50 samples from chr 12
```sh
grep -v "sample_id" /private/groups/migalab/juklucas/censat_regions/active_arrays/asat_arrays_chr12.tsv | shuf -n 50 | cut -f 1-3 | sed 's/\t/,/g' > /private/groups/patenlab/mira/centrolign/batch_submissions/runtime_debugging/r2_QC_v2_flanks_test/chr12/samples.txt
```
2. Regenerate fasta files with 100kb flanks

```sh
# Add 100 kb to bed files
cat /private/groups/patenlab/mira/centrolign/batch_submissions/runtime_debugging/r2_QC_v2_flanks_test/chr12/samples.txt | cut -f1-2 -d"," | sed 's/,/./g' | while read line ; do
  grep chr12 /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/${line}_asat_arrays.bed | awk 'BEGIN{OFS="\t"} {$2=($2-100000<0)?0:$2-100000; $3=$3+100000; print}' > /private/groups/patenlab/mira/centrolign/batch_submissions/runtime_debugging/r2_QC_v2_flanks_test/chr12/per_smp_asat_beds_100kb/${line}_asat_arrays.100kb_flank.bed
done

# Extract fasta files with flanks
git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis pull

cd /private/groups/patenlab/mira/centrolign/batch_submissions/runtime_debugging/r2_QC_v2_flanks_test/chr12/extract_fastas_100kb_flank

mkdir -p logs
sbatch \
  /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/extract_fasta_r2_QCv2_100kb_flank_chr12.sh
```
3. Combine fastas for MSA

```sh
# check none are empty
find . -type f -empty

# combine 100kb fastas
cd /private/groups/patenlab/mira/centrolign/batch_submissions/runtime_debugging/r2_QC_v2_flanks_test/chr12/extract_fastas_100kb_flank

ls | grep fasta | while read line ; do cat $line ; done > /private/groups/patenlab/mira/centrolign/batch_submissions/runtime_debugging/r2_QC_v2_flanks_test/chr12/HPRC_r2_QC_v2_chr12_shuf50_100kb_flanks.fasta

# combine non-flank fastas
cat /private/groups/patenlab/mira/centrolign/batch_submissions/runtime_debugging/r2_QC_v2_flanks_test/chr12/samples.txt | cut -f1-2 -d"," | sed 's/,/./g' | while read line
  do grep $line /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas/fasta_lists/release2_QC_v2_chr12.txt | while read line
    do cat $line
  done
done > /private/groups/patenlab/mira/centrolign/batch_submissions/runtime_debugging/r2_QC_v2_flanks_test/chr12/HPRC_r2_QC_v2_chr12_shuf50.fasta
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/runtime_debugging/r2_QC_v2_flanks_test/chr12/100kb_flanks_MSA/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr12_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/runtime_debugging/r2_QC_v2_flanks_test/chr12/HPRC_r2_QC_v2_chr12_shuf50_100kb_flanks.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/runtime_debugging/r2_QC_v2_flanks_test/chr12/100kb_flanks_MSA/chr12.100kb_flanks.gfa
```

```sh
#!/bin/bash
#SBATCH --job-name=runtime_test_chr12
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
  -S /private/groups/patenlab/mira/centrolign/batch_submissions/runtime_debugging/r2_QC_v2_flanks_test/chr12/asat_MSA/subproblems/ \
  -T /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/nj_trees/chr12_r2_QC_v2_centrolign_all_pairs_nj_tree.nwk \
  /private/groups/patenlab/mira/centrolign/batch_submissions/runtime_debugging/r2_QC_v2_flanks_test/chr12/HPRC_r2_QC_v2_chr12_shuf50.fasta \
  > /private/groups/patenlab/mira/centrolign/batch_submissions/runtime_debugging/r2_QC_v2_flanks_test/chr12/asat_MSA/chr12.gfa
```

```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/runtime_debugging/r2_QC_v2_flanks_test/chr12/asat_MSA
mkdir -p logs
sbatch slurm.sh
```
