### Experiments to test runtime of centrolign

#### Does adding 100kb of flank really increase the runtime drastically?

For chr12, select 50 samples. Run with 100 kb of flanks, and without 100 kb of flanks

1. Randomly select 50 samples from chr 12
```
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

cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/extract_fastas

mkdir -p logs
sbatch \
  /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/extract_fasta_r2_QCv2_100kb_flank_chr12.sh


```
