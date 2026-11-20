#!/bin/bash -l
#$ -l h_rt=6:00:00
#$ -l mem=4G
#$ -pe smp 4
# One task per sample: recal cram -> recal bam (+ index) in $BAM_DIR.
# usage: cram_to_bam.job.sh <list of sample ids>
source "$PIPELINE_DIR/lib/common.sh"
s=$(task_line < "$1")
module load samtools

cram=$(recal_cram "$s")
out=$BAM_DIR/$s.recal.bam
mkdir -p "$BAM_DIR"
log "$s: $cram -> $out"
samtools view -@ 4 -b -T "$REF_FASTA" -o "$out.tmp" "$cram"
mv "$out.tmp" "$out"
samtools index -@ 4 "$out"
