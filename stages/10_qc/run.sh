#!/bin/bash
# Stage 10 -- read and alignment QC, per sample, then one MultiQC report.
# sarek already runs FastQC/mosdepth/samtools stats when started from FASTQs
# (not with --step variant_calling); Qualimap is not part of sarek or mgp1000.
#
#   ./wes.sh qc fastqc     FastQC on each sample's FASTQs (if listed) + recal bam -> QC/fastqc/<S>/
#   ./wes.sh qc qualimap   Qualimap bamqc on each recal bam                       -> QC/qualimap/<S>/
#   ./wes.sh qc multiqc    Qualimap multi-bamqc + MultiQC over QC/ and sarek's reports/
#                                                       -> QC/qualimap/multi_bamqc/, QC/multiqc/
source "$(dirname "$0")/../../lib/common.sh"
here=$PIPELINE_DIR/stages/10_qc

sample_list() {
    local list=$SAMPLESHEET_DIR/qc_samples.txt s
    mkdir -p "$SAMPLESHEET_DIR"
    sample_rows | cut -f2 > "$list"
    while read -r s; do need "$(recal_bam "$s")"; done < "$list"
    echo "$list"
}

step_fastqc() {
    local list; list=$(sample_list)
    submit fastqc "$(wc -l < "$list")" "$here/fastqc.job.sh" "$list"
}

step_qualimap() {
    local list; list=$(sample_list)
    [[ -n $QC_TARGET_BED ]] && need "$QC_TARGET_BED"
    [[ -z $QC_TARGET_BED ]] && warn "QC_TARGET_BED is blank: Qualimap will report whole-genome coverage, not on-target"
    submit qualimap "$(wc -l < "$list")" "$here/qualimap.job.sh" "$list"
}

step_multiqc() {
    compgen -G "$QC_DIR/*/*" > /dev/null || [[ -d $SAREK_DIR/reports ]] ||
        die "nothing to summarise yet -- run qc fastqc / qc qualimap first"
    submit multiqc 1 "$here/multiqc.job.sh"
}

run_steps "usage: wes.sh qc <fastqc|qualimap|multiqc>" "${1:-}" fastqc qualimap multiqc
