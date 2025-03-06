###############################################################################
##                             create input jsons                            ##
###############################################################################
## workflow name = Asm2AsmAlignerPaf

## on personal computer...

# Remove top up data from data table

mkdir -p ~/Desktop/github_repos/centrolign_analysis/batch_submissions/Asm2AsmAlignerPaf/initial_test/Asm2AsmAlignerPaf_input_jsons
cd ~/Desktop/github_repos/centrolign_analysis/batch_submissions/Asm2AsmAlignerPaf/initial_test/Asm2AsmAlignerPaf_input_jsons

python3 /Users/miramastoras/Desktop/Paten_lab/hprc_intermediate_assembly/hpc/launch_from_table.py \
     --data_table ../Asm2AsmAlignerPaf.csv \
     --field_mapping ../Asm2AsmAlignerPaf_input_mapping.csv \
     --workflow_name Asm2AsmAlignerPaf

## add/commit/push to github (hprc_intermediate_assembly)

###############################################################################
##                             create launch workflow                      ##
###############################################################################

## on HPC...

## check that github repo is up to date
git -C  /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/ pull

# move to working dir
mkdir -p /private/groups/patenlab/mira/centrolign/batch_submissions/Asm2AsmAlignerPaf/initial_test/
cd /private/groups/patenlab/mira/centrolign/batch_submissions/Asm2AsmAlignerPaf/initial_test/

## get files
cp -r /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/batch_submissions/Asm2AsmAlignerPaf/initial_test/* ./

mkdir -p slurm_logs
export PYTHONPATH="/private/home/juklucas/miniconda3/envs/toil/bin/python"

# submit job
sbatch \
     --job-name=Asm2AsmAlignerPaf \
     --array=[11-189]%30 \
     --partition=medium \
     --time=12:00:00 \
     --cpus-per-task=32 \
     --exclude=phoenix-[09,10,22,23,24,18] \
     --mem=400gb \
     --mail-type=FAIL,END \
     --mail-user=mmastora@ucsc.edu \
     /private/groups/hprc/hprc_intermediate_assembly/hpc/toil_sbatch_single_machine.sh \
     --wdl ~/progs/hpp_production_workflows/QC/wdl/tasks/Asm2AsmAlignerPaf.wdl \
     --sample_csv Asm2AsmAlignerPaf.csv \
     --input_json_path '../Asm2AsmAlignerPaf_input_jsons/${SAMPLE_ID}_Asm2AsmAlignerPaf.json'

###############################################################################
##                             write output files                   ##
###############################################################################

cd /private/groups/patenlab/mira/centrolign/batch_submissions/Asm2AsmAlignerPaf/initial_test/

python3 /private/groups/hprc/polishing/hprc_intermediate_assembly/hpc/update_table_with_outputs.py \
      --input_data_table Asm2AsmAlignerPaf.csv  \
      --output_data_table Asm2AsmAlignerPaf.results.csv  \
      --json_location '{sample_id}_Asm2AsmAlignerPaf_outputs.json'
