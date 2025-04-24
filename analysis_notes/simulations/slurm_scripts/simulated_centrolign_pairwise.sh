#!/bin/bash
# Slurm script to align simulated sequences pairwise with centrolign, Unialigner, and RaMa
# and to analyze the accuracy, storing the results in a table
#SBATCH --job-name=simulated-centrolign_pairwise
#SBATCH --partition=medium
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --mem=500gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-60
#SBATCH --output=simulated_centrolign_slurm_logs/analysis_array_job_%A_task_%a.log
#SBATCH --time=12:00:00

date
hostname
pwd

# list cases in $SIMDIR/cases.txt
CHR=$1
DATE=20250421

SIMDIR=/private/groups/patenlab/mira/centrolign/simulations/pairwise_simulations/pair_"$CHR"_sim_cases_"$DATE"/
CASE=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$SIMDIR"/cases.txt)
CASEDIR=$SIMDIR/$CASE
WORKDIR=$SIMDIR/work

CENTROLIGN_OUTFILE=$CASEDIR/aln_centrolign.txt
UNIALIGNER_OUTFILE=$CASEDIR/aln_unialigner.txt
RAMA_OUTFILE=$CASEDIR/aln_rama.txt
#WFA_OUTFILE=$CASEDIR/aln_wfa.txt
#WINNOWMAP_OUTFILE=$CASEDIR/aln_winnowmap.txt
#MINIMAP2_OUTFILE=$CASEDIR/aln_minimap2.txt

#WFA=/private/groups/patenlab/jeizenga/GitHub/WFA2-lib/bin/align_benchmark
CENTROLIGN_DIR=/private/home/mmastora/progs/centrolign/build/
SCRIPTS_DIR=/private/groups/patenlab/mira/centrolign/github/centromere-scripts/

# scores: M,X,O1,E1,O2,E2
# heuristic: min distance, max distance from lead, reduction interval
#WFA_PARAMS="-a gap-affine2p-wfa --affine2p-penalties -20,80,100,30,5000,1 --wfa-heuristic wfa-adaptive --wfa-heuristic-parameters 1000,20000,50 --wfa-memory ultralow --wfa-span global"


CENTROLIGN=$CENTROLIGN_DIR/centrolign
TRUTH_COMPARE=$CENTROLIGN_DIR/compare_truth_aln
UNIALIGNER=/private/groups/patenlab/mira/centrolign/github/unialigner/tandem_aligner/build/bin/tandem_aligner
TO_RAW_SEQ=$SCRIPTS_DIR/data_processing_utils/fasta_to_raw_seq.py
ANALYZE_CASE=/private/groups/patenlab/mira/centrolign/github/centrolign_analysis/scripts/analyze_pair_case.py

mkdir -p $WORKDIR
cd $WORKDIR

echo "simulation:" $CASEDIR

FASTA1=${CASEDIR}/sim_seq1.fasta
FASTA2=${CASEDIR}/sim_seq2.fasta

echo "fasta 1:" $FASTA1
echo "fasta 2:" $FASTA2

echo "aligning with centrolign"
TEMP_FASTA=${WORKDIR}/sim_joined_"$SLURM_ARRAY_TASK_ID".fa
mkdir -p `dirname $TEMP_FASTA`
cat $FASTA1 $FASTA2 > $TEMP_FASTA
mkdir -p `dirname $CENTROLIGN_OUTFILE`
ulimit -v 471859200 && /usr/bin/time -v ${CENTROLIGN} -v 3 $TEMP_FASTA > $CENTROLIGN_OUTFILE
#rm $TEMP_FASTA

echo "aligning with unaligner"
mkdir -p `dirname $UNIALIGNER_OUTFILE`
UNIALIGNER_TEMP_OUTDIR=$WORKDIR/uni_tmp_out_"$SLURM_ARRAY_TASK_ID"
ulimit -v 471859200 && /usr/bin/time -v ${UNIALIGNER} --first $FASTA1 --second $FASTA2 -o $UNIALIGNER_TEMP_OUTDIR
mv $UNIALIGNER_TEMP_OUTDIR/cigar.txt $UNIALIGNER_OUTFILE
# delete the rest of the output
#rm -r $UNIALIGNER_TEMP_OUTDIR

echo "aligning with RAMA"
RAMA_TEMP_OUTDIR=$WORKDIR/rama_tmp_out_"$SLURM_ARRAY_TASK_ID"
mkdir -p `dirname $RAMA_OUTFILE`
ulimit -v 471859200 && /usr/bin/time -v \
docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
    miramastoras/rama:latest ./RaMA -t 5 \
    -r $FASTA1 \
    -q $FASTA2 \
    -o $RAMA_TEMP_OUTDIR

mv $RAMA_TEMP_OUTDIR/cigar.txt $RAMA_OUTFILE

# delete the rest of the output
rm -r $RAMA_TEMP_OUTDIR

# echo "aligning with WFA"
# RAW_SEQ_TEMP=${WORKDIR}/tmp_raw_seq_"$SLURM_ARRAY_TASK_ID".txt
# WFA_TEMP_OUT=${WORKDIR}/tmp_wfa_out_"$SLURM_ARRAY_TASK_ID".txt
# $TO_RAW_SEQ $FASTA1 > $RAW_SEQ_TEMP
# $TO_RAW_SEQ $FASTA2 >> $RAW_SEQ_TEMP
# # limit memory to 32 gB and runtime to 30 min, but don't consider it a failure if we don't get it
# true || timeout -v 30m ulimit -m 33554432 $WFA $WFA_PARAMS -i $RAW_SEQ_TEMP -o $WFA_TEMP_OUT
# # remove the score from the output
# if [ -f $WFA_TEMP_OUT ]; then
#    cut -f 2 $WFA_TEMP_OUT > $WFA_OUTFILE
# else
#    touch $WFA_OUTFILE
# fi
# rm -f $RAW_SEQ_TEMP
# rm -f $WFA_TEMP_OUT

# echo "aligning with winnowmap"
#
# WINNOWMAP_TEMP_OUTDIR=$WORKDIR/winnow_tmp_out_"$SLURM_ARRAY_TASK_ID"
#
# mkdir -p ${WINNOWMAP_TEMP_OUTDIR}
# mkdir -p ${WINNOWMAP_OUTFILE}
#
# time docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
#     mobinasri/long_read_aligner:v0.3.3 \
#     meryl count k=19 output ${WINNOWMAP_TEMP_OUTDIR}/merylDB $FASTA1
#
# time docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
#     mobinasri/long_read_aligner:v0.3.3 \
#     meryl print greater-than distinct=0.9998 \
#     ${WINNOWMAP_TEMP_OUTDIR}/merylDB \
#     > ${WINNOWMAP_TEMP_OUTDIR}/repetitive_k19.txt
#
# docker run -u `id -u`:`id -g` -v /private/groups:/private/groups \
#     mobinasri/long_read_aligner:v0.3.3 \
#     winnowmap \
#     -W ${WINNOWMAP_TEMP_OUTDIR}/repetitive_k19.txt \
#     --eqx -cx asm20 -t 1 \
#     $FASTA1 $FASTA2 \
#     > ${WINNOWMAP_OUTFILE}
#
# #rm -r $WINNOWMAP_TEMP_OUTDIR

# do this from outside the directory to get more sensible output
cd $SIMDIR
python3 $ANALYZE_CASE $CASE $TRUTH_COMPARE
