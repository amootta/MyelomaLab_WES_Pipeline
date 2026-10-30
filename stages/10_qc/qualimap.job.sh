#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -l mem=8G
#$ -pe smp 4
# One task per sample: Qualimap bamqc on the recal bam, duplicates skipped.
# From the old qualimap script. With QC_TARGET_BED set, coverage is reported
# inside/outside the capture regions (Qualimap wants a 6-column BED).
# usage: qualimap.job.sh <list of samples>
source "$PIPELINE_DIR/lib/common.sh"
s=$(task_line < "$1")
use_conda qc
unset DISPLAY   # Qualimap is a Java GUI app; keep it headless

target=()
[[ -n $QC_TARGET_BED ]] && target=(-gff "$QC_TARGET_BED")

out=$QC_DIR/qualimap/$s
mkdir -p "$out"
cd "$out"
qualimap bamqc \
    -bam "$(recal_bam "$s")" \
    "${target[@]}" \
    -outdir "$out" \
    -outformat HTML \
    -nt 4 \
    --skip-duplicated \
    --java-mem-size=24G
