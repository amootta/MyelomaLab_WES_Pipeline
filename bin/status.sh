#!/bin/bash
# wes.sh status -- one row per sample, one column per stage output:
# + present, . missing, blank = not applicable (e.g. HLA only on normals).
source "$(dirname "$0")/../lib/common.sh"

mark() { if have "$@"; then printf '%-7s' +; else printf '%-7s' .; fi; }
blank() { printf '%-7s' ""; }

printf '%-26s %-4s %-7s%-7s%-7s%-7s%-7s%-7s%-7s%-7s%-7s%-7s%-7s\n' \
    sample type cram bam mut2f ascat cons.v hlaT1K hlaOpt hlaLA kirT1K igcall pyclon
while IFS=$'\t' read -r p s st _; do
    printf '%-26s %-4s ' "$s" "$([[ $st == 0 ]] && echo N || echo T)"
    mark "$(recal_cram "$s")"
    mark "$(recal_bam "$s")"
    if [[ $st == 1 ]]; then
        n=$(normals | awk -F'\t' -v p="$p" '$1 == p { print $2; exit }')
        if [[ -n $n ]]; then
            mark "$(mutect2_filtered "$s" "$n")"
            mark "$(ascat_pp "$s" "$n")"
        else blank; blank; fi
        if compgen -G "$CONSENSUS_DIR/pvacseq_input/$p/*.vcf" > /dev/null &&
            grep -lq "	$s" "$CONSENSUS_DIR"/pvacseq_input/"$p"/*.vcf 2> /dev/null; then
            printf '%-7s' +
        else printf '%-7s' .; fi
        blank; blank; blank; blank
        mark "$IGCALLER_DIR/$s"
        mark "$PYCLONE_DIR/$p/$p.results.tsv"
    else
        blank; blank; blank
        mark "$HLA_DIR/T1K/${p}_HLA_genotype.tsv"
        if [[ -e $HLA_DIR/calls/optitype/${p}_result.tsv ]] ||
            [[ -n $(find "$HLA_DIR" -path "*optitype*" -name "${p}*_result.tsv" 2> /dev/null | head -1) ]]; then
            printf '%-7s' +
        else printf '%-7s' .; fi
        mark "$HLA_DIR/HLA_LA/$p/hla/R1_bestguess_G.txt"
        mark "$KIR_DIR/T1K/${p}_KIR_genotype.tsv"
        blank; blank
    fi
    echo
done < <(sample_rows)

echo
echo "pVACseq: $(ls -d "$PVACSEQ_DIR"/*/ 2> /dev/null | grep -vc '/hla/$' || true) patient run folder(s) in $PVACSEQ_DIR"
echo "Logs:    $LOG_DIR"
