#!/bin/bash
# wes.sh -- single entry point for the Myeloma Lab WES pipeline.
#
#   ./wes.sh check                     are the config paths and tools in place?
#   ./wes.sh status                    which outputs exist, per sample
#   ./wes.sh <stage> <step>            submit a step (see docs/STAGES.md)
#
#   DRY_RUN=1 ./wes.sh <stage> <step>  print the qsub line instead of submitting
#   PIPELINE_CONFIG=other.config ./wes.sh ...   use another batch's config
set -eo pipefail
here=$(cd "$(dirname "$0")" && pwd)

usage() {
    cat >&2 << 'EOF'
usage: ./wes.sh <command> [step]

  check                                     config paths + tools
  status                                    per-sample progress table

  01 sarek       fastq | cram | bams | filter_mutect2
  02 mgp1000     somatic
  03 consensus   vcf | isec
  04 hla         t1k | hlatyping | hlala | collect | concordance
  05 kir         t1k | kircle | concordance
  06 igcaller    run
  07 pvacseq     hla | run
  08 pyclone     run
  09 signatures  input

Steps are submitted one at a time: each checks the previous step's
output exists before queueing anything. Order: docs/STAGES.md
EOF
    exit 1
}

cmd=${1:-}
shift || true
case $cmd in
    check)  exec bash "$here/bin/check.sh" "$@" ;;
    status) exec bash "$here/bin/status.sh" "$@" ;;
    sarek|01)      stage=01_sarek ;;
    mgp1000|02)    stage=02_mgp1000 ;;
    consensus|03)  stage=03_consensus ;;
    hla|04)        stage=04_hla ;;
    kir|05)        stage=05_kir ;;
    igcaller|06)   stage=06_igcaller ;;
    pvacseq|07)    stage=07_pvacseq ;;
    pyclone|08)    stage=08_pyclone ;;
    signatures|09) stage=09_signatures ;;
    *) usage ;;
esac
exec bash "$here/stages/$stage/run.sh" "$@"
