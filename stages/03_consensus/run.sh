#!/bin/bash
# Stage 03 -- turn mgp1000's three-caller output into the files later stages use.
#
#   ./wes.sh consensus vcf     union consensus TSVs -> PASS in >= N callers -> one VCF
#                              per tumour, sorted into pvacseq_input/<patient>/ (stage 07)
#   ./wes.sh consensus isec    per-caller VCFs -> normalise -> hard filter -> sites in
#                              >= 2 of 3 callers -> variant tables for MMSig (stage 09)
source "$(dirname "$0")/../../lib/common.sh"
here=$PIPELINE_DIR/stages/03_consensus

# mgp1000's output is shared by every batch; both steps only take pairs whose
# tumour is in this batch's samples.tsv.
step_vcf() {
    local f pair n=0 t
    compgen -G "$MGP1000_CONSENSUS_GLOB" > /dev/null || die "no files match MGP1000_CONSENSUS_GLOB=$MGP1000_CONSENSUS_GLOB -- has stage 02 finished?"
    for f in $MGP1000_CONSENSUS_GLOB; do
        pair=$(basename "$f")
        t=$(pair_tumour "${pair%%.hq.union.*}")
        [[ -n $t ]] && n=$((n + 1))
    done
    ((n > 0)) || die "none of the mgp1000 consensus tables belong to tumours in $SAMPLES"
    log "$n consensus table(s) belong to this batch"
    while IFS=$'\t' read -r _ t; do
        compgen -G "$MGP1000_SOMATIC_DIR/consensus/${t}[_.]*" > /dev/null || warn "$t: no mgp1000 consensus table yet"
    done < <(tumours)
    submit consensus_vcf 1 "$here/consensus_vcf.job.sh"
}

step_isec() {
    local c
    for c in mutect strelka varscan; do
        [[ -d $MGP1000_SOMATIC_DIR/$c ]] || die "missing $MGP1000_SOMATIC_DIR/$c -- has stage 02 finished?"
    done
    submit consensus_isec 1 "$here/isec.job.sh"
}

run_steps "usage: wes.sh consensus <vcf|isec>" "${1:-}" vcf isec
