#!/bin/bash
#SBATCH --job-name=msa_simulations_direct_pairwise
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=56gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=[1-510]%128
#SBATCH --output=logs/array_job_%A_task_%a.log
#SBATCH --time=1:00:00


/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/active_array_csvs/test.csv

samtools faidx -r $REGIONFILE ~{assemblyFasta} | sed "s/>/>$SAMPLE.$PARNUM /g" > $HORFASTA
