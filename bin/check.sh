#!/bin/bash
# wes.sh check -- confirm every path in the config exists and the sample sheet
# is well formed, before anything is submitted.
source "$(dirname "$0")/../lib/common.sh"

bad=0
check() {  # check <label> <path>
    if [[ -e $2 ]]; then printf '  ok       %-22s %s\n' "$1" "$2"
    else printf '  MISSING  %-22s %s\n' "$1" "$2"; bad=1; fi
}

echo "config:  $PIPELINE_CONFIG"
echo "samples: $SAMPLES"
echo
echo "reference + resources"
for v in REF_FASTA REF_FAI DBSNP KNOWN_INDELS GNOMAD_AF PON FUNCOTATOR_DS ASCAT_LOCI ASCAT_ALLELES \
         ASCAT_LOCI_GC ASCAT_LOCI_RT MAPPABILITY CONTROLFREEC_CHR_DIR; do
    check "$v" "${!v}"
done
echo
echo "tools"
for v in CONDA_SH T1K_DIR KIRCLE_DIR HLATYPING_DIR HLAHD_TARBALL IGCALLER_HOME IGCALLER_SAMTOOLS MGP1000_DIR \
         PVACSEQ_PIPELINE PVACSEQ_IEDB PVACSEQ_VEP_PLUGINS; do
    check "$v" "${!v}"
done
command -v nextflow > /dev/null && echo "  ok       nextflow               $(command -v nextflow)" \
    || { echo "  MISSING  nextflow               (not on PATH)"; bad=1; }
command -v qsub > /dev/null || echo "  note     qsub not on PATH -- run this on a Myriad login node"

echo
echo "sample sheet"
awk -F'\t' '
    NR == 1 { n = NF; for (i = 1; i <= NF; i++) c[$i] = i
              split("patient sample status sex lane fastq_1 fastq_2 cram bam", want, " ")
              for (w in want) if (!(want[w] in c)) { print "  MISSING column: " want[w]; bad = 1 }
              next }
    NF != n { print "  line " NR ": " NF " columns, header has " n; bad = 1 }
    $c["status"] != 0 && $c["status"] != 1 { print "  line " NR ": status must be 0 (normal) or 1 (tumour)"; bad = 1 }
    $c["sex"] != "XX" && $c["sex"] != "XY" { print "  line " NR ": sex must be XX or XY"; bad = 1 }
    $c["fastq_1"] == "" && $c["cram"] == "" && $c["bam"] == "" { print "  line " NR ": no fastq, cram or bam"; bad = 1 }
    END { exit bad }
' "$SAMPLES" || bad=1
printf '  %d patients, %d normals, %d tumours, %d tumour/normal pairs\n' \
    "$(patients | wc -l)" "$(normals | wc -l)" "$(tumours | wc -l)" "$(pairs 2> /dev/null | wc -l)"
pairs > /dev/null   # print any "no normal" warnings

echo
if ((bad)); then echo "Problems found -- fix config/pipeline.config or the sample sheet."; exit 1; fi
echo "All good."
