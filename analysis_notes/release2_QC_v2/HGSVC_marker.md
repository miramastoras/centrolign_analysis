### Code for comparing connectedness curves for HGSVC and HPRC r2 for the marker 

#### Run all pairs centrolign alignments for HGSVC 

Concatenate all of the QC csv files together. Have to run chr3 and chr4 separately because they came in later.

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
  /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/extract_fasta_HPRC_marker.sh
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
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/centrolign_direct_all_pairs.sh \
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
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/centrolign_direct_all_pairs.sh \
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
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/centrolign_direct_all_pairs.sh \
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
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/centrolign_direct_all_pairs.sh \
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
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/centrolign_direct_all_pairs.sh \
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
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/centrolign_direct_all_pairs.sh \
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
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/centrolign_direct_all_pairs.sh \
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
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/centrolign_direct_all_pairs.sh \
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
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/centrolign_direct_all_pairs.sh \
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
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/centrolign_direct_all_pairs.sh \
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
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/centrolign_direct_all_pairs.sh \
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
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/centrolign_direct_all_pairs.sh \
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
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/centrolign_direct_all_pairs.sh \
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
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/centrolign_direct_all_pairs.sh \
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
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/centrolign_direct_all_pairs.sh \
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
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/centrolign_direct_all_pairs.sh \
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
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/centrolign_direct_all_pairs.sh \
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
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/centrolign_direct_all_pairs.sh \
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
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/centrolign_direct_all_pairs.sh \
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
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/centrolign_direct_all_pairs.sh \
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
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/centrolign_direct_all_pairs.sh \
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
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/centrolign_direct_all_pairs.sh \
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
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/centrolign_direct_all_pairs.sh \
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
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/centrolign_direct_all_pairs.sh \
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
    /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/centrolign_direct_all_pairs.sh \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}/hgsvc_all_pairs_combinations_${CHR}.txt \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/hgsvc_all_pairs/${CHR}

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
  /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/analysis_notes/release2_QC_v2/slurm_scripts/extract_fasta_HPRC_smps_marker.sh
```