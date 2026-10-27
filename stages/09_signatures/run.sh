#!/bin/bash
# Stage 09 -- mutational signatures with MMSig (R, run interactively).
#
#   ./wes.sh signatures input   merge stage 03's per-pair 2-of-3 tables into
#                               consensus/isec/mmsig_input.tsv (sample chr pos ref alt)
#
# Then in R (https://github.com/evenrus/mmsig):
#   muts <- read.delim("mmsig_input.tsv")
#   ... mm_fit_signatures(muts.input = muts, ...)
# NOTE: these calls are hg38, and mmsig assumes hg19 by default -- set the
#       genome option when you run it.
source "$(dirname "$0")/../../lib/common.sh"

step_input() {
    local in=$CONSENSUS_DIR/isec/variant_tsv out=$CONSENSUS_DIR/isec/mmsig_input.tsv
    compgen -G "$in/*.tsv" > /dev/null || die "no variant tables in $in -- run 'wes.sh consensus isec' first"
    {
        echo -e "sample\tchr\tpos\tref\talt"
        for f in "$in"/*.tsv; do tail -n +2 "$f" | cut -f1-5; done
    } > "$out"
    log "$out: $(($(wc -l < "$out") - 1)) SNVs from $(ls "$in"/*.tsv | wc -l) pair(s)"
}

run_steps "usage: wes.sh signatures input" "${1:-}" input
