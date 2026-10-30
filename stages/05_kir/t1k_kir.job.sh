#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -l mem=4G
#$ -pe smp 4
# One task per patient: T1K KIR genotyping on the normal bam.
# usage: t1k_kir.job.sh <normals.tsv: patient normal>
source "$PIPELINE_DIR/lib/common.sh"
IFS=$'\t' read -r p n < <(task_line < "$1")
use_conda t1k

mkdir -p "$KIR_DIR/T1K"
run-t1k -b "$(recal_bam "$n")" \
    --preset kir-wes \
    -f "$T1K_DIR/kiridx/kiridx_dna_seq.fa" \
    -c "$T1K_DIR/kiridx/kiridx_dna_coord.fa" \
    -o "$KIR_DIR/T1K/${p}_KIR" \
    -t 4
