### Code for comparing connectedness curves for HGSVC and HPRC r2 for the marker 

#### Run all pairs centrolign alignments for HGSVC 

Concatenate all of the QC csv files together. 

```sh
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

grep project /private/groups/hprc/qc_testing/hgsvc/censat_region_comparison/arrays_by_chromosome/chr1/asat-hgsvc-pass.csv > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_asat_arrays.csv

for chr in "${chromosomes[@]}"
do
  wc -l /private/groups/hprc/qc_testing/hgsvc/censat_region_comparison/arrays_by_chromosome/${chr}/asat-hgsvc-pass.csv
  cat /private/groups/hprc/qc_testing/hgsvc/censat_region_comparison/arrays_by_chromosome/${chr}/asat-hgsvc-pass.csv | grep -v "project" >> /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_asat_arrays.csv
done  
```

Python script to create per sample bed files containing all of the arrays
```sh
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/parse_QC_csv_marker.py \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_asat_arrays.csv \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_per_smp_asat_beds/
```

Get list of all samples for submitting sbatch script
```sh
grep -v project /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_asat_arrays.csv | cut -f2-3 -d","  | sed 's/\t/,/g' > hgsvc_all_samples_with_asats.txt

sort hgsvc_all_samples_with_asats.txt | uniq > tmp ; mv tmp hgsvc_all_samples_with_asats.txt
```

Slurm script to run on list of samples, downloads fasta file per sample and extracts HOR sequence per sample, placing it in dir per chromosome
```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/extract_fastas_hgsvc/

mkdir -p logs

sbatch \
  /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/extract_fasta_HPRC_marker.sh
```

Accounting check: Make sure we aren't missing any fastas
```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/extract_fastas_hgsvc/

# Check for empty fastas
echo "=== Empty fastas ==="
find . -type f -empty
echo "(none above = good)"

# Compare unique samples in CSV vs fastas per chrom
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

echo ""
echo "chrom,fastas,unique_samples_in_csv,match,missing"

for chr in "${chromosomes[@]}"; do
    fasta_dir="/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/extract_fastas_hgsvc/${chr}"
    csv="/private/groups/hprc/qc_testing/hgsvc/censat_region_comparison/arrays_by_chromosome/${chr}/asat-hgsvc-pass.csv"

    fastas=$(ls "${fasta_dir}/" 2>/dev/null | wc -l)

    # Unique sample_haplotype (col 4) in CSV
    unique_csv=$(tail -n +2 "${csv}" | cut -d',' -f4 | sort -u | wc -l)

    match=$([[ "$fastas" -eq "$unique_csv" ]] && echo "TRUE" || echo "MISMATCH")

    # Find which samples are missing
    missing=$(comm -23 \
        <(tail -n +2 "${csv}" | cut -d',' -f4 | sort -u) \
        <(ls "${fasta_dir}/" | sed "s/_${chr}_hor_array.fasta//" | sort) \
        | paste -sd ',' -)

    echo "${chr},${fastas},${unique_csv},${match},${missing}"
done

```

Get list of fasta files per chromosome
```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/extract_fastas_hgsvc/

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/extract_fastas_hgsvc/fasta_lists/

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")
for chr in "${chromosomes[@]}"
do
    ls /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/extract_fastas_hgsvc/${chr}/ | while read line ;
      do realpath $chr/$line >> /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/extract_fastas_hgsvc/fasta_lists/hgsvc_${chr}.txt
    done
done
```

Generate all vs all combinations
```sh
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

# get all vs all pairwise combinations
for chr in "${chromosomes[@]}"
do
  mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${chr}

  while read -r s1; do
      while read -r s2; do
          [[ "$s1" < "$s2" ]] && echo -e "$s1\t$s2\t$chr"
      done < /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/extract_fastas_hgsvc/fasta_lists/hgsvc_${chr}.txt
  done < /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/extract_fastas_hgsvc/fasta_lists/hgsvc_${chr}.txt > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${chr}/hgsvc_all_pairs_combinations_${chr}.txt
done

for chr in "${chromosomes[@]}"
do
  wc -l /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${chr}/hgsvc_all_pairs_combinations_${chr}.txt
done
```

Run centrolign all pairs 
```sh
CHR="chr1"
ARRAY="[1-7875]%128"

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/

sbatch \
    --job-name=${CHR}_hgsvc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/hgsvc_all_pairs_combinations_${CHR}.txt \ # combinations file
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR} # chrom dir
```

Claude generated commands:
```sh
CHR="chr1"
ARRAY="1-7875%128"

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/

sbatch \
    --job-name=${CHR}_hgsvc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/hgsvc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}

CHR="chr2"
ARRAY="1-6786%128"

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/

sbatch \
    --job-name=${CHR}_hgsvc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/hgsvc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}

CHR="chr3"
ARRAY="1-7021%128"

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/

sbatch \
    --job-name=${CHR}_hgsvc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/hgsvc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}

CHR="chr4"
ARRAY="1-5886%128"

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/

sbatch \
    --job-name=${CHR}_hgsvc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/hgsvc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}

CHR="chr5"
ARRAY="1-7503%128"

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/

sbatch \
    --job-name=${CHR}_hgsvc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/hgsvc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}

CHR="chr6"
ARRAY="1-6105%128"

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/

sbatch \
    --job-name=${CHR}_hgsvc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/hgsvc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}

CHR="chr7"
ARRAY="1-7503%128"

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/

sbatch \
    --job-name=${CHR}_hgsvc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/hgsvc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}

#### Stopped here 

CHR="chr8"
ARRAY="1-7626%128"

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/

sbatch \
    --job-name=${CHR}_hgsvc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/hgsvc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}

CHR="chr9"
ARRAY="1-7626%128"

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/

sbatch \
    --job-name=${CHR}_hgsvc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/hgsvc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}

CHR="chr10"
ARRAY="1-7875%128"

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/

sbatch \
    --job-name=${CHR}_hgsvc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/hgsvc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}

CHR="chr11"
ARRAY="1-8001%128"

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/

sbatch \
    --job-name=${CHR}_hgsvc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/hgsvc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}

CHR="chr12"
ARRAY="1-7503%128"

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/

sbatch \
    --job-name=${CHR}_hgsvc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/hgsvc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}

CHR="chr13"
ARRAY="1-7626%128"

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/

sbatch \
    --job-name=${CHR}_hgsvc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/hgsvc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}

CHR="chr14"
ARRAY="1-7503%128"

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/

sbatch \
    --job-name=${CHR}_hgsvc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/hgsvc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}

CHR="chr15"
ARRAY="1-7626%128"

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/

sbatch \
    --job-name=${CHR}_hgsvc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/hgsvc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}

CHR="chr16"
ARRAY="1-7626%128"

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/

sbatch \
    --job-name=${CHR}_hgsvc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/hgsvc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}

CHR="chr17"
ARRAY="1-6670%128"

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/

sbatch \
    --job-name=${CHR}_hgsvc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/hgsvc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}

CHR="chr18"
ARRAY="1-6903%128"

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/

sbatch \
    --job-name=${CHR}_hgsvc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/hgsvc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}

CHR="chr19"
ARRAY="1-7875%128"

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/

sbatch \
    --job-name=${CHR}_hgsvc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/hgsvc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}

CHR="chr20"
ARRAY="1-7140%128"

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/

sbatch \
    --job-name=${CHR}_hgsvc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/hgsvc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}

CHR="chr21"
ARRAY="1-7140%128"

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/

sbatch \
    --job-name=${CHR}_hgsvc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/hgsvc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}

CHR="chr22"
ARRAY="1-7750%128"

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/

sbatch \
    --job-name=${CHR}_hgsvc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/hgsvc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}

CHR="chrX"
ARRAY="1-4186%128"

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/

sbatch \
    --job-name=${CHR}_hgsvc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/hgsvc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}

CHR="chrY"
ARRAY="1-351%128"

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/

sbatch \
    --job-name=${CHR}_hgsvc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/hgsvc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}

```
Create distance matrices 

Get distance matrices for all pairs
```sh
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/distance_matrices

for chr in "${chromosomes[@]}"
do
    python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/cigar_to_distance.py -d 1 \
    -a /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${chr}/pairwise_cigar/ \
    -o /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/distance_matrices/${chr}_hgsvc_centrolign_direct
  done
```

#### Run all pairs centrolign alignments for HPRC  

Concatenate all of the QC csv files together. Have to run chr3 and chr4 separately because they came in later.

```sh
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

grep project /private/groups/hprc/qc_testing/hgsvc/censat_region_comparison/arrays_by_chromosome/chr1/asat-hprc-pass.csv > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_asat_arrays.csv

for chr in "${chromosomes[@]}"
do
  wc -l /private/groups/hprc/qc_testing/hgsvc/censat_region_comparison/arrays_by_chromosome/${chr}/asat-hprc-pass.csv
  cat /private/groups/hprc/qc_testing/hgsvc/censat_region_comparison/arrays_by_chromosome/${chr}/asat-hprc-pass.csv | grep -v "project" >> /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_asat_arrays.csv
done  
```

Python script to create per sample bed files containing all of the arrays
```sh
python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/parse_QC_csv_marker.py \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_asat_arrays.csv \
  /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_per_smp_asat_beds/
```

Get list of all samples for submitting sbatch script
```sh
cat /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_samples_with_asats.txt | grep -v "HG06807,1" | grep -v "HG06807,2" | grep -v "CHM13,0" | grep -v "HG002,1" | grep -v "HG002,2" > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_samples_with_asats.txt
```

Slurm script to run on list of samples, downloads fasta file per sample and extracts HOR sequence per sample, placing it in dir per chromosome
```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/extract_fastas_hprc

mkdir -p logs

sbatch \
  /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/extract_fasta_HPRC_smps_marker.sh
```

Accounting check: Make sure we aren't missing any fastas
```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/extract_fastas_hprc/

# Check for empty fastas
echo "=== Empty fastas ==="
find . -type f -empty
echo "(none above = good)"

# Compare unique samples in CSV vs fastas per chrom
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

echo ""
echo "chrom,fastas,unique_samples_in_csv,match,missing"

for chr in "${chromosomes[@]}"; do
    fasta_dir="/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/extract_fastas_hprc/${chr}"
    csv="/private/groups/hprc/qc_testing/hgsvc/censat_region_comparison/arrays_by_chromosome/${chr}/asat-hprc-pass.csv"

    fastas=$(ls "${fasta_dir}/" 2>/dev/null | wc -l)

    # Unique sample_haplotype (col 4) in CSV
    unique_csv=$(tail -n +2 "${csv}" | cut -d',' -f4 | sort -u | wc -l)

    match=$([[ "$fastas" -eq "$unique_csv" ]] && echo "TRUE" || echo "MISMATCH")

    # Find which samples are missing
    missing=$(comm -23 \
        <(tail -n +2 "${csv}" | cut -d',' -f4 | sort -u) \
        <(ls "${fasta_dir}/" | sed "s/_${chr}_hor_array.fasta//" | sort) \
        | paste -sd ',' -)

    echo "${chr},${fastas},${unique_csv},${match},${missing}"
done

```

Get list of fasta files per chromosome
```sh
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/extract_fastas_hprc/

mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/extract_fastas_hprc/fasta_lists/

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")
for chr in "${chromosomes[@]}"
do
    ls /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/extract_fastas_hprc/${chr}/ | while read line ;
      do realpath $chr/$line >> /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/extract_fastas_hprc/fasta_lists/hprc_${chr}.txt
    done
done
```
Generate all vs all combinations
```sh
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

# get all vs all pairwise combinations
for chr in "${chromosomes[@]}"
do
  mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${chr}

  while read -r s1; do
      while read -r s2; do
          [[ "$s1" < "$s2" ]] && echo -e "$s1\t$s2\t$chr"
      done < /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/extract_fastas_hprc/fasta_lists/hprc_${chr}.txt
  done < /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/extract_fastas_hprc/fasta_lists/hprc_${chr}.txt > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${chr}/hprc_all_pairs_combinations_${chr}.txt
done

for chr in "${chromosomes[@]}"
do
  wc -l /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${chr}/hprc_all_pairs_combinations_${chr}.txt
done
```

Run all pairs 
```sh
# chr1 (99,681 pairs)
CHR="chr1"
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/

ARRAY="1-30000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="30001-60000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="60001-90000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="90001-99681%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

# chr2 (91,806 pairs)
CHR="chr2"
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/

ARRAY="1-30000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="30001-60000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="60001-90000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="90001-91806%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

# chr3 (57,291 pairs)
CHR="chr3"
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/

ARRAY="1-30000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="30001-57291%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

# chr4 (61,776 pairs)
CHR="chr4"
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/

ARRAY="1-30000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="30001-60000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="60001-61776%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

# chr5 (100,576 pairs)
CHR="chr5"
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/

ARRAY="1-30000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="30001-60000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="60001-90000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="90001-100576%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

# chr6 (62,481 pairs)
CHR="chr6"
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/

ARRAY="1-30000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="30001-60000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="60001-62481%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

# chr7 (91,378 pairs)
CHR="chr7"
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/

ARRAY="1-30000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="30001-60000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="60001-90000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="90001-91378%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

# chr8 (96,580 pairs)
CHR="chr8"
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/

ARRAY="1-30000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="30001-60000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="60001-90000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="90001-96580%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

# chr9 (99,235 pairs)
CHR="chr9"
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/

ARRAY="1-30000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="30001-60000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="60001-90000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="90001-99235%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

# chr10 (102,378 pairs)
CHR="chr10"
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/

ARRAY="1-30000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="30001-60000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="60001-90000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="90001-102378%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

# chr11 (102,831 pairs)
CHR="chr11"
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/

ARRAY="1-30000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="30001-60000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="60001-90000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="90001-102831%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

# chr12 (101,475 pairs)
CHR="chr12"
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/

ARRAY="1-30000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="30001-60000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="60001-90000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="90001-101475%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

# chr13 (78,210 pairs)
CHR="chr13"
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/

ARRAY="1-30000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="30001-60000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="60001-78210%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

# chr14 (96,580 pairs)
CHR="chr14"
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/

ARRAY="1-30000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="30001-60000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="60001-90000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="90001-96580%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

# chr15 (86,320 pairs)
CHR="chr15"
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/

ARRAY="1-30000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="30001-60000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="60001-86320%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

# chr16 (97,903 pairs)
CHR="chr16"
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/

ARRAY="1-30000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="30001-60000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="60001-90000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="90001-97903%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

# chr17 (39,060 pairs)
CHR="chr17"
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/

ARRAY="1-30000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="30001-39060%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

# chr18 (34,980 pairs)
CHR="chr18"
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/

ARRAY="1-30000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="30001-34980%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

# chr19 (100,128 pairs)
CHR="chr19"
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/

ARRAY="1-30000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="30001-60000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="60001-90000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="90001-100128%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

# chr20 (94,830 pairs)
CHR="chr20"
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/

ARRAY="1-30000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="30001-60000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="60001-90000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="90001-94830%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

# chr21 (75,466 pairs)
CHR="chr21"
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/

ARRAY="1-30000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="30001-60000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="60001-75466%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

# chr22 (94,395 pairs)
CHR="chr22"
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/

ARRAY="1-30000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="30001-60000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="60001-90000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="90001-94395%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

# chrX (55,945 pairs)
CHR="chrX"
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/

ARRAY="1-30000%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

ARRAY="30001-55945%128"
sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

# chrY (2,016 pairs)
CHR="chrY"
ARRAY="1-2016%128"
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/logs/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/

sbatch \
    --job-name=${CHR}_hprc_all_pairs \
    --array=$ARRAY \
    --export=CHR=${CHR} \
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/HPRC_marker/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}/hprc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hprc_all_pairs/${CHR}

```

#### Pairwise tree heatmap with CenHap tree 

chr 11 
```sh
/private/groups/migalab/juklucas/centrolign/cenhap_assignment/cenhap_inference_out/chr1/

docker run -u `id -u`:`id -g` \
    -v /private/groups:/private/groups/ \
    miramastoras/tree_heatmap:latest \
    python3 /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/pairwise_tree_heatmap_v2.py \
    -t /private/groups/migalab/juklucas/centrolign/cenhap_assignment/cenhap_inference_out/chr7/chr7.mc_common_sites.upgma.nwk \
    -s /private/groups/patenlab/mira/chr7_samples.txt \
    -p /private/groups/patenlab/mira/chr7_dist_matrix.csv \
    -m "Centrolign pairwise match distance" \
    -n "chr7 CenHap tree" \
    -d "Pairwise Match Distance" \
    -o /private/groups/patenlab/mira/chr7_cenhap_tree_ \
    --no_labels
```

Pairwise tree heatmap:

Tree and sample lists from: /private/groups/migalab/juklucas/centrolign/triangle_heatmaps/2026_03_12_r2_triangle_heatmaps/cenhap_tree/chr11

```
chr=chr11
python3 /Users/miramastoras/Desktop/github_repos/centrolign_analysis/scripts/pairwise_tree_heatmap_v2.py \
    -t /Users/miramastoras/Desktop/pruned_tree.nwk \
    -s /Users/miramastoras/Desktop/chr11_samples.txt \
    -p /Users/miramastoras/Desktop/HPRC_release2_QCv2_all_pairs_heatmaps/${chr}_r2_QC_v2_centrolign_pairwise_distance.csv \
    -m "Centrolign pairwise match distance" \
    -n "" \
    -d "CenHap Tree" \
    -o /Users/miramastoras/Desktop/${chr}_r2_QC_v2_cenhap2 --no_labels \
    --cenhap_labels /Users/miramastoras/Desktop/chr11.cenhap_predictions.tsv \
    --cenhap_colors /Users/miramastoras/Desktop/cenhap_colors.tsv
```

