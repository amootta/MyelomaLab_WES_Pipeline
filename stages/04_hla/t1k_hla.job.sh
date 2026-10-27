#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -l mem=4G
#$ -pe smp 4
# One task per patient: T1K HLA genotyping on the normal bam.
# Fixes the old t1k_hla.sh, which pointed at the KIR index (kiridx) and wrote
# into KIR/ -- so its "HLA" output was really a second KIR run.
# usage: t1k_hla.job.sh <normals.tsv: patient normal>
source "$PIPELINE_DIR/lib/common.sh"
IFS=$'\t' read -r p n < <(task_line < "$1")
use_conda t1k

mkdir -p "$HLA_DIR/T1K"
# TODO: T1K's presets are hla / hla-wgs / kir-wgs / kir-wes. hla-wgs is what was
#       used before; for exome data plain `hla` may be the better fit.
run-t1k -b "$(recal_bam "$n")" \
    --preset hla-wgs \
    -f "$T1K_DIR/hlaidx/hlaidx_dna_seq.fa" \
    -c "$T1K_DIR/hlaidx/hlaidx_dna_coord.fa" \
    -o "$HLA_DIR/T1K/${p}_HLA" \
    -t 4
