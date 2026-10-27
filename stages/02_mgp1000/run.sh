#!/bin/bash
# Stage 02 -- MGP1000 somatic calling (Mutect2 + Strelka2 + VarScan consensus)
# on the recal bams from stage 01. mgp1000 only reads bams, so run
# `wes.sh sarek bams` first.
#
#   ./wes.sh mgp1000 somatic
source "$(dirname "$0")/../../lib/common.sh"
here=$PIPELINE_DIR/stages/02_mgp1000

# TODO: mgp1000_WES.sh is referenced in the old README but was never
#       committed. Check that the input folder below is where your mgp1000
#       checkout reads bams from, and paste your real command into
#       mgp1000.job.sh.
MGP1000_BAM_DIR=$MGP1000_DIR/input/preprocessedBams

step_somatic() {
    local sheet=$SAMPLESHEET_DIR/mgp1000.csv p t n tb nb
    mkdir -p "$SAMPLESHEET_DIR" "$MGP1000_BAM_DIR"
    echo "normal,tumor" > "$sheet"
    while IFS=$'\t' read -r p t n; do
        tb=$(recal_bam "$t") nb=$(recal_bam "$n")
        need "$tb" "$nb"
        # mgp1000 takes bare filenames resolved against its own input folder
        ln -sf "$tb" "$MGP1000_BAM_DIR/$t.bam"
        ln -sf "$nb" "$MGP1000_BAM_DIR/$n.bam"
        for b in "$tb:$t" "$nb:$n"; do
            local src=${b%%:*} name=${b##*:}
            if [[ -e $src.bai ]]; then ln -sf "$src.bai" "$MGP1000_BAM_DIR/$name.bam.bai"
            elif [[ -e ${src%.bam}.bai ]]; then ln -sf "${src%.bam}.bai" "$MGP1000_BAM_DIR/$name.bam.bai"
            else die "no index for $src"; fi
        done
        echo "$n.bam,$t.bam" >> "$sheet"
    done < <(pairs)
    log "wrote $sheet; bams linked into $MGP1000_BAM_DIR"
    submit mgp1000_somatic 1 "$here/mgp1000.job.sh" "$sheet"
}

run_steps "usage: wes.sh mgp1000 somatic" "${1:-}" somatic
