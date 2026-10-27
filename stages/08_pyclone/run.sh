#!/bin/bash
# Stage 08 -- PyClone-VI clonal structure across each patient's tumour timepoints.
# Inputs: mgp1000 2-of-3 consensus SNV VCFs (stage 03 vcf), ASCAT cnvs +
# purity/ploidy (stage 01 sarek), recal crams (stage 01).
#
#   ./wes.sh pyclone run    one job per patient, all steps:
#        union of consensus SNVs across timepoints (build_pyclone_input.R --emit-sites)
#        -> forced allele counts in every sample (alleleCounter)
#        -> PyClone-VI input table               (build_pyclone_input.R)
#        -> pyclone-vi fit + write-results-file
#        -> cluster table + CCF plot             (summarize_pyclone.R)
#   -> pyclone/<patient>/
source "$(dirname "$0")/../../lib/common.sh"
here=$PIPELINE_DIR/stages/08_pyclone

step_run() {
    local list=$SAMPLESHEET_DIR/pyclone_patients.txt p t n dir man vcf cnv pp
    mkdir -p "$SAMPLESHEET_DIR"
    : > "$list"
    while read -r p; do
        dir=$PYCLONE_DIR/$p
        man=$dir/manifest.tsv
        mkdir -p "$dir"
        echo -e "CASE_ID\tSAMPLE\tpass_vcf\tcram\tcounts\tascat_cnv\tascat_pp\tsex" > "$man"
        while IFS=$'\t' read -r _ t n; do
            vcf=$(consensus_vcf "$t") cnv=$(ascat_cnv "$t" "$n") pp=$(ascat_pp "$t" "$n")
            [[ -n $vcf ]] || die "$t: no consensus VCF in $CONSENSUS_DIR/vcf -- run 'wes.sh consensus vcf' first"
            need "$vcf" "$cnv" "$pp" "$(recal_cram "$t")"
            echo -e "$p\t$t\t$vcf\t$(recal_cram "$t")\t$dir/counts/$t.forced.tsv\t$cnv\t$pp\t$(sex_of "$p")" >> "$man"
        done < <(pairs | awk -F'\t' -v p="$p" '$1 == p')
        if (($(wc -l < "$man") > 1)); then echo "$p" >> "$list"; else rm -rf "$dir"; fi
    done < <(patients)
    log "manifests written under $PYCLONE_DIR/<patient>/manifest.tsv"
    submit pyclone "$(wc -l < "$list")" "$here/pyclone.job.sh" "$list"
}

run_steps "usage: wes.sh pyclone run" "${1:-}" run
