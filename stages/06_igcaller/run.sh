#!/bin/bash
# Stage 06 -- IgCaller: IG rearrangements, class switch and translocations,
# per tumour/normal pair.
#
#   ./wes.sh igcaller run      -> IgCaller/<tumour>/
source "$(dirname "$0")/../../lib/common.sh"
here=$PIPELINE_DIR/stages/06_igcaller

step_run() {
    local list=$SAMPLESHEET_DIR/pairs.tsv p t n
    mkdir -p "$SAMPLESHEET_DIR"
    pairs > "$list"
    while IFS=$'\t' read -r p t n; do need "$(mapped_bam "$t")" "$(mapped_bam "$n")"; done < "$list"
    need "$IGCALLER_HOME/IgCaller" "$IGCALLER_HOME/IgCaller_reference_files"
    submit igcaller "$(wc -l < "$list")" "$here/igcaller.job.sh" "$list"
}

run_steps "usage: wes.sh igcaller run" "${1:-}" run
