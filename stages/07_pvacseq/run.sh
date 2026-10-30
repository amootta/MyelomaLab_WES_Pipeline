#!/bin/bash
# Stage 07 -- neoantigen prediction with nextflow-pvacseq, one run per patient.
# Joins stage 03 (consensus VCFs) with stage 04 (HLA typing).
#
#   ./wes.sh pvacseq hla    write pvacseq/hla/<patient>.csv from the HLA calls
#                           (caller chosen by HLA_SOURCE_FOR_PVACSEQ) -- check it!
#   ./wes.sh pvacseq run    one nextflow-pvacseq head job per patient -> pvacseq/<patient>/
source "$(dirname "$0")/../../lib/common.sh"
here=$PIPELINE_DIR/stages/07_pvacseq

step_hla() {
    need "$CONSENSUS_DIR/pvacseq_input" "$HLA_DIR/calls/$HLA_SOURCE_FOR_PVACSEQ"
    module load "$PYTHON_MODULE" 2> /dev/null || true
    python3 "$PIPELINE_DIR/bin/build_pvacseq_hla.py" \
        --calls-dir "$HLA_DIR/calls" \
        --source "$HLA_SOURCE_FOR_PVACSEQ" \
        --vcf-dir "$CONSENSUS_DIR/pvacseq_input" \
        --outdir "$PVACSEQ_DIR/hla"
}

step_run() {
    local list=$SAMPLESHEET_DIR/pvacseq_patients.txt p
    mkdir -p "$SAMPLESHEET_DIR"
    : > "$list"
    for d in "$CONSENSUS_DIR"/pvacseq_input/*/; do
        p=$(basename "$d")
        [[ -e $PVACSEQ_DIR/hla/$p.csv ]] || { warn "$p: no HLA csv, run 'wes.sh pvacseq hla' -- skipped"; continue; }
        echo "$p" >> "$list"
    done
    need "$PVACSEQ_PIPELINE" "$PVACSEQ_IEDB" "$PVACSEQ_VEP_PLUGINS"
    submit pvacseq "$(wc -l < "$list")" "$here/pvacseq.job.sh" "$list"
}

run_steps "usage: wes.sh pvacseq <hla|run>" "${1:-}" hla run
