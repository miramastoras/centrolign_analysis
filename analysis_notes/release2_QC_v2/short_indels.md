### Calling short indels from pairwise cigar strings

Run short indel calling script
```sh
cd /private/groups/patenlab/mira/centrolign/analysis/short_indels_pairwise

mkdir -p logs

sbatch --array=1 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/call_short_indels_pairwise.sh \
    /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/12012025_completed_subgroups.csv \
    /private/groups/patenlab/mira/centrolign/analysis/short_indels_pairwise
```
