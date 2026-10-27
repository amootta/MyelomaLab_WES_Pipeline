#!/bin/bash -l
#$ -l h_rt=5:00:00
#$ -l mem=40G
#$ -pe smp 1
# One task per tumour/normal pair: re-filter sarek's raw Mutect2 calls with
# the lab's thresholds, then annotate PASS calls with Funcotator (MAF).
# Merges the old filter_annotate.sh and mutect2_filter_annotate.sh.
# usage: filter_mutect2.job.sh <pairs.tsv: patient tumour normal>
source "$PIPELINE_DIR/lib/common.sh"
IFS=$'\t' read -r _ t n < <(task_line < "$1")
module load gatk

id=$(pair_id "$t" "$n")
raw=$(mutect2_raw "$t" "$n")
out=$(mutect2_filtered "$t" "$n")
dir=$(dirname "$raw")
mkdir -p "$MUTECT2_DIR"

# Use sarek's contamination / segmentation / orientation-bias outputs when it
# wrote them -- FilterMutectCalls is weaker without them, and
# build_pyclone_input.R warns when --ob-priors was skipped.
extra=()
[[ -e $dir/$id.mutect2.contamination.table ]] && extra+=(--contamination-table "$dir/$id.mutect2.contamination.table")
[[ -e $dir/$id.mutect2.segmentation.table ]]  && extra+=(--tumor-segmentation "$dir/$id.mutect2.segmentation.table")
[[ -e $dir/$id.mutect2.artifactprior.tar.gz ]] && extra+=(--ob-priors "$dir/$id.mutect2.artifactprior.tar.gz")
log "$id: extra filter inputs: ${extra[*]:-none}"
# TODO: sarek already writes $id.mutect2.filtered.vcf.gz with its own
#       FilterMutectCalls run. This re-filter only adds the 0.01 VAF floor and
#       the MQ/BQ 20 thresholds -- confirm that is still wanted.

gatk FilterMutectCalls \
    -R "$REF_FASTA" \
    -V "$raw" \
    -O "$out" \
    --min-median-mapping-quality 20 \
    --min-median-base-quality 20 \
    --min-allele-fraction 0.01 \
    "${extra[@]}"

gatk Funcotator \
    -R "$REF_FASTA" \
    -V "$out" \
    -O "${out%.vcf.gz}.maf" \
    --remove-filtered-variants true \
    --three-prime-flank-size 1000 \
    --output-file-format MAF \
    --ref-version hg38 \
    --annotation-default "tumor_barcode:$t" \
    --data-sources-path "$FUNCOTATOR_DS"
