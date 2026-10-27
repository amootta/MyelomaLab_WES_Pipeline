#!/bin/bash
# Stage 02 -- MGP1000 somatic calling (Mutect2 + Strelka2 + VarScan consensus,
# plus Manta/Delly/IgCaller) on the recal bams from stage 01. mgp1000 only
# reads bams, so run `wes.sh sarek bams` first.
#
#   ./wes.sh mgp1000 somatic
#
# Starting from FASTQs instead (mgp1000's own preprocessing.nf, not sarek's
# bams): see the commented preprocessing.nf command in mgp1000.job.sh.
source "$(dirname "$0")/../../lib/common.sh"
here=$PIPELINE_DIR/stages/02_mgp1000

# Put every pair's bam + index into one --input_dir. Hard links, not symlinks:
# mgp1000's containers only see the bound input folder, so a symlink pointing
# elsewhere on Scratch would dangle inside them. Falls back to a copy across
# filesystems.
stage_bam() {  # stage_bam <sample>
    local src dst idx
    src=$(recal_bam "$1") dst=$MGP1000_INPUT_DIR/$1.bam
    need "$src"
    [[ $src -ef $dst ]] || ln -f "$src" "$dst" 2> /dev/null || cp "$src" "$dst"
    for idx in "$src.bai" "${src%.bam}.bai"; do
        if [[ -e $idx ]]; then
            [[ $idx -ef $dst.bai ]] || ln -f "$idx" "$dst.bai" 2> /dev/null || cp "$idx" "$dst.bai"
            return
        fi
    done
    die "no index for $src"
}

step_somatic() {
    local sheet=$SAMPLESHEET_DIR/mgp1000_somatic.csv t n
    mkdir -p "$SAMPLESHEET_DIR" "$MGP1000_INPUT_DIR"
    echo "normal,tumor" > "$sheet"
    while IFS=$'\t' read -r _ t n; do
        stage_bam "$t"
        stage_bam "$n"
        echo "$n.bam,$t.bam" >> "$sheet"
    done < <(pairs)
    log "wrote $sheet; bams in $MGP1000_INPUT_DIR"
    submit mgp1000_somatic 1 "$here/mgp1000.job.sh" "$sheet"
}

run_steps "usage: wes.sh mgp1000 somatic" "${1:-}" somatic
