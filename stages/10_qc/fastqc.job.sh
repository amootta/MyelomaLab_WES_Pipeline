#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -l mem=4G
#$ -pe smp 2
# One task per sample: FastQC on its FASTQs (every lane listed in samples.tsv
# that still exists) and on its recal bam.
# usage: fastqc.job.sh <list of samples>
source "$PIPELINE_DIR/lib/common.sh"
s=$(task_line < "$1")
use_conda qc

files=()
while read -r f; do
    if [[ -e $f ]]; then files+=("$f"); else warn "$s: fastq $f not found, skipped"; fi
done < <(awk -F'\t' -v s="$s" '
    NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
    $c["sample"] == s { if ($c["fastq_1"] != "") print $c["fastq_1"]; if ($c["fastq_2"] != "") print $c["fastq_2"] }
' "$SAMPLES")
files+=("$(recal_bam "$s")")

out=$QC_DIR/fastqc/$s
mkdir -p "$out"
log "$s: FastQC on ${#files[@]} file(s)"
fastqc -t 2 -o "$out" "${files[@]}"
