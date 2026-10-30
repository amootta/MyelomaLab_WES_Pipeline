#!/bin/bash -l
#$ -l h_rt=40:00:00
#$ -l mem=16G
#$ -pe smp 1
# One task per patient: KIRCLE on the normal bam.
# usage: kircle.job.sh <normals.tsv: patient normal>
source "$PIPELINE_DIR/lib/common.sh"
IFS=$'\t' read -r p n < <(task_line < "$1")
module unload gcc-libs
module load blast+
module load python3
module load samtools
export BLASTDB=$KIRCLE_DIR/ref_files/ref_fastas/

# Run from the batch's KIR folder so any scratch files land there, not in the
# shared install.
mkdir -p "$KIR_DIR/KIRCLE"
cd "$KIR_DIR/KIRCLE"
python "$KIRCLE_DIR/KIRCLE/KIRCLE.py" \
    -i "$(recal_bam "$n")" \
    -o "$KIR_DIR/KIRCLE/$p" \
    -l "$KIRCLE_DIR/ref_files/hg38_KIR_locations_noKV_noalt.tsv" \
    -g hg38
