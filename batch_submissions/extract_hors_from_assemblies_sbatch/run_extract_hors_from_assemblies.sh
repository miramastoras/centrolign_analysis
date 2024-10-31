# on hpc
cd /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_from_assemblies_sbatch

git -C  /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/ pull

cp -r /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/batch_submissions/extract_hors_from_assemblies_sbatch/* ./

mkdir -p slurm_logs

sbatch extract_hors_from_assemblies.sh extract_hors_from_assemblies_sbatch.csv 
