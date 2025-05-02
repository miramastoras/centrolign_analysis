#!/bin/bash
# Slurm script to use centrolign as a pairwise aligner on all pairs of sequences
# from an input directory

#SBATCH --job-name=CHM13_pairwise-centrolign
#SBATCH --partition=short
#SBATCH --mail-user=mmastora@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=56gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/array_job_%x_%j_%A_%a.log
#SBATCH --time=1:00:00

set -e

# Fetch input arguments with this while loop
while [ $# -gt 0 ]; do
  case "$1" in
    -h|--help)
     cat << 'EOF'

      Usage:

      sbatch centrolign_all_pairs_CHM13.sh \
        --chr chrY \
        --job-name=chrY_all_pairs_CHM13 \
        --array=[23-26]%4

EOF
      exit 0
      ;;
    -c|--chr)
      shift
      if [ $# -gt 0 ]; then
        export CHR="$1"
      else
        echo "Error: No CHR specified"
        exit 1
      fi
      shift
      ;;
    *)
      break
      ;;
  esac
done

set -x

date
hostname
pwd

CHROMDIR=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2/all_pairs/${CHR}
FASTADIR=/private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/
WORKDIR=$CHROMDIR/work/
OUTDIR=$CHROMDIR/pairwise_cigar/

mkdir -p $OUTDIR
mkdir -p $WORKDIR
mkdir -p $FASTADIR

cd $WORKDIR

SAMPLE1=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$CHROMDIR"/HPRC_release2_contiguous_HOR_CHM13_combinations_${CHR}.txt | cut -f1)
SAMPLE2=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$CHROMDIR"/HPRC_release2_contiguous_HOR_CHM13_combinations_${CHR}.txt | cut -f2)
echo "sample 1:" $SAMPLE1
echo "sample 2:" $SAMPLE2
echo "out:" $OUTDIR

FASTA1="/private/groups/patenlab/mira/centrolign/annotations/chm13/per_chrom_for_centrolign/chm13v2.0.${CHR}.active_hor.upper.fa"
FASTA2=`grep ${SAMPLE2} /private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/HPRC_release2_contiguous_HORs_${CHR}.fasta_list.txt`
FASTA2="/private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/${FASTA2}"

echo "fasta 1:" $FASTA1
echo "fasta 2:" $FASTA2
TEMP_FASTA=${WORKDIR}/${SAMPLE1}_${SAMPLE2}.fa
cat $FASTA1 $FASTA2 > $TEMP_FASTA
time /private/home/mmastora/progs/centrolign/build/centrolign -v 3 --skip-calibration $TEMP_FASTA > $OUTDIR/pairwise_cigar_${SAMPLE1}_${SAMPLE2}.txt
rm $TEMP_FASTA
