#!/bin/bash -l
#$ -l h_rt=4:00:00
#$ -l mem=8G
#$ -pe smp 1
# Qualimap multi-bamqc across every sample that has a bamqc result, then one
# MultiQC report over QC/ and sarek's reports/ (mosdepth, samtools stats, ...).
# sarek's own FastQC is left out: stage 10 reruns FastQC on the same FASTQs,
# and the duplicate sample names would overwrite each other in the report.
# usage: multiqc.job.sh
source "$PIPELINE_DIR/lib/common.sh"
use_conda qc
unset DISPLAY

qm=$QC_DIR/qualimap
input=$qm/multibamqc_input.txt
: > "$input"
for d in "$qm"/*/; do
    d=${d%/}
    [[ -e $d/genome_results.txt ]] && printf '%s\t%s\n' "$(basename "$d")" "$d" >> "$input"
done
if [[ -s $input ]]; then
    log "Qualimap multi-bamqc over $(wc -l < "$input") sample(s)"
    qualimap multi-bamqc -d "$input" -outdir "$qm/multi_bamqc" -outformat HTML --java-mem-size=6G
else
    warn "no Qualimap bamqc results in $qm, skipping multi-bamqc"
fi

dirs=("$QC_DIR")
[[ -d $SAREK_DIR/reports ]] && dirs+=("$SAREK_DIR/reports")
cd "$QC_DIR"
# extra_fn_clean_exts: "<S>.recal" -> "<S>", so the bam's FastQC, Qualimap and
# samtools rows line up under one sample name.
multiqc "${dirs[@]}" \
    --ignore "*/reports/fastqc/*" \
    --cl-config "extra_fn_clean_exts: ['.recal']" \
    --title "$BATCH_NAME" \
    -o "$QC_DIR/multiqc" \
    -n multiqc_report.html \
    -f
