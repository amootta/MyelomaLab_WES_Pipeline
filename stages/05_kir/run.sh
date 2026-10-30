#!/bin/bash
# Stage 05 -- germline KIR typing on each patient's NORMAL bam.
#
#   ./wes.sh kir t1k           T1K (kiridx, kir-wes)  -> KIR/T1K/<patient>_KIR_genotype.tsv
#   ./wes.sh kir kircle        KIRCLE                 -> KIR/KIRCLE/<patient>/
#   ./wes.sh kir concordance   KIR/KIR_concordance.xlsx (T1K vs KIR_TRUTH, or the
#                              6 samples embedded in bin/kir_concordance.py)
source "$(dirname "$0")/../../lib/common.sh"
here=$PIPELINE_DIR/stages/05_kir

normal_list() {
    local list=$SAMPLESHEET_DIR/normals.tsv p n
    mkdir -p "$SAMPLESHEET_DIR"
    normals > "$list"
    while IFS=$'\t' read -r p n; do need "$(recal_bam "$n")"; done < "$list"
    echo "$list"
}

step_t1k() {
    local list; list=$(normal_list)
    need "$T1K_DIR/kiridx/kiridx_dna_seq.fa" "$T1K_DIR/kiridx/kiridx_dna_coord.fa"
    submit t1k_kir "$(wc -l < "$list")" "$here/t1k_kir.job.sh" "$list"
}

step_kircle() {
    local list; list=$(normal_list)
    need "$KIRCLE_DIR/KIRCLE/KIRCLE.py" "$KIRCLE_DIR/ref_files/hg38_KIR_locations_noKV_noalt.tsv"
    submit kircle "$(wc -l < "$list")" "$here/kircle.job.sh" "$list"
}

step_concordance() {
    compgen -G "$KIR_DIR/T1K/*_KIR_genotype.tsv" > /dev/null || die "no T1K KIR results in $KIR_DIR/T1K -- run 'wes.sh kir t1k' first"
    [[ -z $KIR_TRUTH ]] || need "$KIR_TRUTH"
    submit kir_concordance 1 "$here/concordance.job.sh"
}

run_steps "usage: wes.sh kir <t1k|kircle|concordance>" "${1:-}" t1k kircle concordance
