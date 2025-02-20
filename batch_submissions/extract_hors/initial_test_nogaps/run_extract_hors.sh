###############################################################################
##                             create input jsons                            ##
###############################################################################
## workflow name = extract_hors
## experiment = initial_test_nogaps

## on personal computer...

# Remove top up data from data table

mkdir -p ~/Desktop/github_repos/centrolign_analysis/batch_submissions/extract_hors/initial_test_nogaps/extract_hors_input_jsons
cd ~/Desktop/github_repos/centrolign_analysis/batch_submissions/extract_hors/initial_test_nogaps/extract_hors_input_jsons

python3 /Users/miramastoras/Desktop/Paten_lab/hprc_intermediate_assembly/hpc/launch_from_table.py \
     --data_table ../extract_hors_initial_test_nogaps.csv \
     --field_mapping ../extract_hors_input_mapping.csv \
     --workflow_name extract_hors

## add/commit/push to github (hprc_intermediate_assembly)

###############################################################################
##                             create launch workflow                      ##
###############################################################################

## on HPC...

## check that github repo is up to date
git -C  /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/ pull

# move to working dir
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors/initial_test_nogaps
cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors/initial_test_nogaps

## get files
cp -r /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/batch_submissions/extract_hors/initial_test_nogaps/* ./

mkdir -p slurm_logs
export PYTHONPATH="/private/home/juklucas/miniconda3/envs/toil/bin/python"

# submit job
sbatch \
     --job-name=extract_hors \
     --array=[137]%30 \
     --partition=short \
     --time=1:00:00 \
     --cpus-per-task=32 \
     --exclude=phoenix-[09,10,22,23,24,18] \
     --mem=400gb \
     --mail-type=FAIL,END \
     --mail-user=mmastora@ucsc.edu \
     ~/progs/hprc_intermediate_assembly/hpc/toil_sbatch_single_machine_debug.sh \
     --wdl /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/wdl/tasks/extract_hors.wdl \
     --sample_csv extract_hors_initial_test_nogaps.csv \
     --input_json_path '../extract_hors_input_jsons/${SAMPLE_ID}_extract_hors.json'

###############################################################################
##                             write output files                   ##
###############################################################################

cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors/initial_test_nogaps

python3 /private/groups/hprc/polishing/hprc_intermediate_assembly/hpc/update_table_with_outputs.py \
      --input_data_table extract_hors_initial_test_nogaps.csv  \
      --output_data_table extract_hors_initial_test_nogaps.results.csv  \
      --json_location '{sample_id}_extract_hors_outputs.json'
