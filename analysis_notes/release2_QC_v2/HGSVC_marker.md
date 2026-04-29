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
# check to make sure none of the fastas are empty
cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/extract_fastas_hgsvc/

find . -type f -empty

# check that number of fastas per chrom matches julian's original files
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

for chr in "${chromosomes[@]}"
do
  fastas=`ls /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/HPRC_marker/extract_fastas_hgsvc/${chr}/ | wc -l`
  QC=`grep -v "project" /private/groups/hprc/qc_testing/hgsvc/censat_region_comparison/arrays_by_chromosome/${chr}/asat-hgsvc-pass.csv | sort | uniq | wc -l`

  echo $chr,$fastas,$QC
  if [[ "$fastas" == "$QC" ]]; then
    echo "true"
  fi
done

## checks passed.
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