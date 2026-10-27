#!/bin/bash -l
#$ -l h_rt=40:00:00
#$ -l mem=16G
#$ -pe smp 4
# One task per patient: HLA-LA on the normal bam.
# Needs graphs/PRG_MHC_GRCh38_withIMGT copied from the RDSS into the conda
# env (<env>/opt/hla-la/graphs/) -- see the README install notes.
# usage: hla_la.job.sh <normals.tsv: patient normal>
source "$PIPELINE_DIR/lib/common.sh"
IFS=$'\t' read -r p n < <(task_line < "$1")
use_conda hla-la

mkdir -p "$HLA_DIR/HLA_LA"
HLA-LA.pl --BAM "$(recal_bam "$n")" \
    --graph PRG_MHC_GRCh38_withIMGT \
    --sampleID "$p" \
    --maxThreads 4 \
    --workingDir "$HLA_DIR/HLA_LA"
