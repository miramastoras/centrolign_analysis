### Run pairwise tree heatmap on all samples for chr12

Get pairwise cigars
```
#!/bin/bash
#SBATCH --job-name=centrolign_chr12_all_samples_pairwise
#SBATCH --partition=long
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=500gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=centrolign.log
#SBATCH --time=72:00:00


/private/home/mmastora/progs/centrolign/build/centrolign -v 4 \
    -S /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/chr12/jobstore \
    -R \
    -T /private/groups/patenlab/jeizenga/centromere/chr12/KGP4_TRIOS_MAC5_chr12_CPR_EHet30_no_PS_PID_PGT_lifted_over.v1.1_mask.nwk.txt \
    -A /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/chr12/pairwise_cigars/ \
    /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/chr12/chr12_hprc_r2_initial_test_inside_tree.fasta \
    > /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/initial_test/chr12/initial_test_chr12.centrolign.gfa
```
