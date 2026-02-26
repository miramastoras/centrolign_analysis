- for every pair:
- Take CDR coordinates for each sample
- calculate mutation rate in 1,000 bp windows inside CDR
- randomly shuffle 10,000 windows? across rest of the array, calculate mutation rate
- compare the distribution


mkdir permuted_beds


for i in {1..10000}; do
  bedtools shuffle \
    -i cdr_windows.bed \
    -g genome.sizes \
    -incl centromere_minus_CDR.bed \
    -noOverlapping \
    > permuted_beds/perm_${i}.bed
done

### Looking at some synteny plots

```sh
python /private/groups/migalab/juklucas/centrolign/chr12_test125/synteny_plot_bokeh.py   \
    --beds \
        /private/groups/migalab/juklucas/censat_paper/local_identity/2025_11_13/chrX/run/HG03516_2_CM089243.1_57905064_59678080/analysis/local_identity_outputs/e2241006-0cde-4f9c-8223-713f8f691d4a/HG03516_2_CM089243_1_57905064_59678080_local_identity_subset_coords.bed \
        /private/groups/migalab/juklucas/censat_paper/local_identity/2025_11_13/chrX/run/NA18505_2_JBHIJZ010000010.1_5810287_7563244/analysis/local_identity_outputs/548b6c0c-04f2-4157-89b1-f6d9f8b5a7f9/NA18505_2_JBHIJZ010000010_1_5810287_7563244_local_identity_subset_coords.bed \
        /private/groups/migalab/juklucas/censat_paper/local_identity/2025_11_13/chrX/run/NA19700_2_CM089646.1_57846226_59633981/analysis/local_identity_outputs/ed7e6060-3800-4f16-aca8-f95731c4752d/NA19700_2_CM089646_1_57846226_59633981_local_identity_subset_coords.bed \
        /private/groups/migalab/juklucas/censat_paper/local_identity/2025_11_13/chrX/run/HG03579_2_JAGYVT020000063.1_57947258_59737060/analysis/local_identity_outputs/6b55272a-0522-45c6-87fd-126af3540e95/HG03579_2_JAGYVT020000063_1_57947258_59737060_local_identity_subset_coords.bed \
    --cigars \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrX/induced_pairwise_cigars/pairwise_cigar_HG03516.2_NA18505.2.txt \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrX/induced_pairwise_cigars/pairwise_cigar_NA18505.2_NA19700.2.txt \
        /private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chrX/induced_pairwise_cigars/pairwise_cigar_NA19700.2_HG03579.2.txt \
    --output /private/groups/patenlab/mira/chrX_synteny.html \
    --show-mismatches \
    --web
```
