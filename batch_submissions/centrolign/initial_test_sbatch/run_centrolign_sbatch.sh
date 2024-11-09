cd /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test

mkdir -p sbatch_logs

git -C /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/ pull

cp -r /private/groups/patenlab/mira/centrolign/github/centrolign_analysis/batch_submissions/centrolign/initial_test_sbatch/* .

sbatch centrolign_sbatch.sh \
    centrolign_initial_test_sbatch.csv
