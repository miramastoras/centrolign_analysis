### Calling short indels from pairwise cigar strings

Run short indel calling script
```sh
cd /private/groups/patenlab/mira/centrolign/analysis/short_indels_pairwise

mkdir -p logs

sbatch --array=[1-43]%43 \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/call_short_indels_pairwise.sh \
    /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/12012025_completed_subgroups.csv \
    /private/groups/patenlab/mira/centrolign/analysis/short_indels_pairwise
```

Create csv file formatted as clade,chr,path_to_short_indel_beds
```sh
cd /private/groups/patenlab/mira/centrolign/analysis/short_indels_pairwise

echo "clade,chr,path_to_SV_beds" > 12142025_clade_chr_short_indel_beds.csv
cat /private/groups/patenlab/mira/centrolign/analysis/SVs_pairwise/12012025_completed_subgroups.csv | while IFS=',' read -r clade cigars sample_lists ; do
  CHR=$(echo $clade | cut -f1 -d"_")
  echo $clade,$CHR,/private/groups/patenlab/mira/centrolign/analysis/short_indels_pairwise/${CHR}/short_indel_beds/${clade}/ >>  12142025_clade_chr_short_indel_beds.csv
done
```
Get short indels in ASM coords 
