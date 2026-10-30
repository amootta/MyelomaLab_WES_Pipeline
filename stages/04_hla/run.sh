#!/bin/bash
# Stage 04 -- germline HLA typing on each patient's NORMAL bam, three ways,
# then gather the calls in one place for the concordance report and pVACseq.
#
#   ./wes.sh hla t1k          T1K (hlaidx)                      -> HLA/T1K/<patient>_HLA_genotype.tsv
#   ./wes.sh hla hlatyping    patched nf-core/hlatyping (OptiType, HLA-HD) -> HLA/hlatyping/
#   ./wes.sh hla hlala        HLA-LA                            -> HLA/HLA_LA/<patient>/
#   ./wes.sh hla collect      link every caller's per-patient file into HLA/calls/<caller>/
#   ./wes.sh hla concordance  HLA/HLA_typing_comparison.xlsx (needs HLA_TRUTH in the config)
#
# Results are keyed by PATIENT (one normal each), which is how stage 07 looks them up.
source "$(dirname "$0")/../../lib/common.sh"
here=$PIPELINE_DIR/stages/04_hla

normal_list() {
    local list=$SAMPLESHEET_DIR/normals.tsv p n
    mkdir -p "$SAMPLESHEET_DIR"
    normals > "$list"
    while IFS=$'\t' read -r p n; do need "$(recal_bam "$n")"; done < "$list"
    echo "$list"
}

step_t1k() {
    local list; list=$(normal_list)
    need "$T1K_DIR/hlaidx/hlaidx_dna_seq.fa" "$T1K_DIR/hlaidx/hlaidx_dna_coord.fa"
    submit t1k_hla "$(wc -l < "$list")" "$here/t1k_hla.job.sh" "$list"
}

step_hlatyping() {
    local list sheet=$SAMPLESHEET_DIR/hlatyping.csv p n
    list=$(normal_list)
    need "$HLATYPING_DIR/main.nf" "$HLAHD_TARBALL"
    echo "sample,bam,seq_type" > "$sheet"
    while IFS=$'\t' read -r p n; do echo "$p,$(recal_bam "$n"),dna"; done < "$list" >> "$sheet"
    log "wrote $sheet"
    submit hlatyping 1 "$here/hlatyping.job.sh" "$sheet"
}

step_hlala() {
    local list; list=$(normal_list)
    submit hla_la "$(wc -l < "$list")" "$here/hla_la.job.sh" "$list"
}

# One flat folder per caller, one file per patient, named so that
# hla_concordance.py's default sample-ID regexes match.
step_collect() {
    local calls=$HLA_DIR/calls f p n=0
    mkdir -p "$calls"/{optitype,hlahd,t1k,hlala}
    while read -r p; do
        f=$HLA_DIR/T1K/${p}_HLA_genotype.tsv
        [[ -e $f ]] && ln -sf "$f" "$calls/t1k/${p}_HLA_genotype.tsv" && n=$((n + 1))
        # all of HLA_DIR, not just hlatyping/: earlier runs used --outdir HLA/
        f=$(find "$HLA_DIR" -path "$calls" -prune -o -path "*optitype*" -name "${p}*_result.tsv" -print 2>/dev/null | head -1)
        [[ -n $f ]] && ln -sf "$f" "$calls/optitype/${p}_result.tsv" && n=$((n + 1))
        # HLA-HD writes <sample>/result/<sample>_final.result.txt
        f=$(find "$HLA_DIR" -path "$calls" -prune -o -path "*hlahd*" -name "${p}*final*.txt" -print 2>/dev/null | head -1)
        [[ -n $f ]] && ln -sf "$f" "$calls/hlahd/${p}_hlahd_final.txt" && n=$((n + 1))
        f=$HLA_DIR/HLA_LA/$p/hla/R1_bestguess_G.txt
        [[ -e $f ]] && ln -sf "$f" "$calls/hlala/${p}_bestguess_G.txt" && n=$((n + 1))
    done < <(normals | cut -f1)
    log "linked $n call files into $calls"
    for f in "$calls"/*; do echo "  $(basename "$f"): $(ls "$f" | wc -l) patient(s)"; done
}

step_concordance() {
    [[ -n $HLA_TRUTH ]] || die "set HLA_TRUTH in the config (bin/hla_concordance.py --truth-template makes a blank one)"
    need "$HLA_TRUTH" "$HLA_DIR/calls"
    submit hla_concordance 1 "$here/concordance.job.sh"
}

run_steps "usage: wes.sh hla <t1k|hlatyping|hlala|collect|concordance>" "${1:-}" t1k hlatyping hlala collect concordance
