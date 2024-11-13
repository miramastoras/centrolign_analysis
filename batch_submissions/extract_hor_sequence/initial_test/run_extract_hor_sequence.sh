###############################################################################
##                             create input jsons                            ##
###############################################################################
## workflow name = extract_hor_sequence
## experiment = initial_test

## on personal computer...

# Remove top up data from data table

mkdir -p ~/Desktop/github_repos/centrolign_analysis/batch_submissions/extract_hor_sequence/initial_test/extract_hor_sequence_input_jsons
cd ~/Desktop/github_repos/centrolign_analysis/batch_submissions/extract_hor_sequence/initial_test/extract_hor_sequence_input_jsons

python3 /Users/miramastoras/Desktop/Paten_lab/hprc_intermediate_assembly/hpc/launch_from_table.py \
     --data_table ../extract_hor_sequence_initial_test.csv \
     --field_mapping ../extract_hor_sequence_input_mapping.csv \
     --workflow_name extract_hor_sequence

## add/commit/push to github (hprc_intermediate_assembly)

###############################################################################
##                             create launch workflow                      ##
###############################################################################

## on HPC...

## check that github repo is up to date
git -C  /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/ pull

# move to working dir
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hor_sequence/initial_test/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hor_sequence/initial_test/

## get files
cp -r /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/batch_submissions/extract_hor_sequence/initial_test/* ./

mkdir -p slurm_logs
export PYTHONPATH="/private/home/juklucas/miniconda3/envs/toil/bin/python"

# submit job
sbatch \
     --job-name=extract_hor_sequence \
     --array=[1-189]%30 \
     --partition=short \
     --time=1:00:00 \
     --cpus-per-task=32 \
     --exclude=phoenix-[09,10,22,23,24,18] \
     --mem=400gb \
     --mail-type=FAIL,END \
     --mail-user=mmastora@ucsc.edu \
     /private/groups/hprc/hprc_intermediate_assembly/hpc/toil_sbatch_single_machine.sh \
     --wdl /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/wdl/tasks/extract_hor_sequence.wdl \
     --sample_csv extract_hor_sequence_initial_test.csv \
     --input_json_path '../extract_hor_sequence_input_jsons/${SAMPLE_ID}_extract_hor_sequence.json'

###############################################################################
##                             write output files                   ##
###############################################################################

cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hor_sequence/initial_test

python3 /private/groups/hprc/polishing/hprc_intermediate_assembly/hpc/update_table_with_outputs.py \
      --input_data_table extract_hor_sequence_initial_test.csv  \
      --output_data_table extract_hor_sequence_initial_test.results.csv  \
      --json_location '{sample_id}_extract_hor_sequence_outputs.json'
