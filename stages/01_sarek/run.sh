#!/bin/bash
# Stage 01 -- nf-core/sarek: alignment, BQSR, Mutect2/Strelka/Manta/ASCAT/CNVkit,
# then the extra Mutect2 filtering + Funcotator annotation.
#
#   ./wes.sh sarek fastq          full run from FASTQs (mapping onwards)
#   ./wes.sh sarek cram           variant calling only, from existing recal crams
#   ./wes.sh sarek bams           cram -> bam for every sample (HLA/KIR/IgCaller/mgp1000 need bams)
#   ./wes.sh sarek filter_mutect2 FilterMutectCalls + Funcotator MAF per tumour/normal pair
source "$(dirname "$0")/../../lib/common.sh"
here=$PIPELINE_DIR/stages/01_sarek

# Every row of a sample is one lane: sarek maps each lane as its own read
# group, then merges them into one bam per sample before MarkDuplicates/BQSR.
# Read groups are named sample + lane, so a lane number repeated within a
# sample (re-sequenced on another flowcell) gets a _2, _3 ... suffix.
step_fastq() {
    local sheet=$SAMPLESHEET_DIR/sarek_fastq.csv
    mkdir -p "$SAMPLESHEET_DIR"
    awk -F'\t' 'BEGIN { OFS = "," }
        NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i
                  print "patient", "sex", "status", "sample", "lane", "fastq_1", "fastq_2"; next }
        $c["fastq_1"] != "" {
            lane = $c["lane"]
            if (++seen[$c["sample"] SUBSEP lane] > 1) {
                lane = lane "_" seen[$c["sample"] SUBSEP lane]
                print "NOTE: " $c["sample"] " lane " $c["lane"] " listed more than once; this row becomes lane " lane > "/dev/stderr"
            }
            print $c["patient"], $c["sex"], $c["status"], $c["sample"], lane, $c["fastq_1"], $c["fastq_2"]
        }
    ' "$SAMPLES" > "$sheet"
    need $(tail -n +2 "$sheet" | cut -d, -f6,7 | tr , ' ')
    log "wrote $sheet ($(($(wc -l < "$sheet") - 1)) lanes)"
    submit sarek_fastq 1 "$here/sarek.job.sh" fastq "$sheet"
}

step_cram() {
    local sheet=$SAMPLESHEET_DIR/sarek_cram.csv s p st sx cram
    mkdir -p "$SAMPLESHEET_DIR"
    echo "patient,sex,status,sample,cram,crai" > "$sheet"
    while IFS=$'\t' read -r p s st sx; do
        cram=$(recal_cram "$s")
        [[ -n $(sample_field "$s" bam) && -z $(sample_field "$s" cram) ]] && continue  # bam-only samples
        need "$cram" "$cram.crai"
        echo "$p,$sx,$st,$s,$cram,$cram.crai" >> "$sheet"
    done < <(sample_rows)
    log "wrote $sheet"
    submit sarek_cram 1 "$here/sarek.job.sh" cram "$sheet"
}

step_bams() {
    local list=$SAMPLESHEET_DIR/cram_to_bam.txt s
    mkdir -p "$SAMPLESHEET_DIR"
    : > "$list"
    while IFS=$'\t' read -r _ s _ _; do
        [[ -e $(recal_bam "$s") ]] && continue   # override, sarek's own, or already converted
        need "$(recal_cram "$s")"
        echo "$s" >> "$list"
    done < <(sample_rows)
    submit cram_to_bam "$(wc -l < "$list")" "$here/cram_to_bam.job.sh" "$list"
}

step_filter_mutect2() {
    local list=$SAMPLESHEET_DIR/pairs.tsv p t n
    mkdir -p "$SAMPLESHEET_DIR"
    pairs > "$list"
    while IFS=$'\t' read -r p t n; do need "$(mutect2_raw "$t" "$n")"; done < "$list"
    submit filter_mutect2 "$(wc -l < "$list")" "$here/filter_mutect2.job.sh" "$list"
}

run_steps "usage: wes.sh sarek <fastq|cram|bams|filter_mutect2>" "${1:-}" fastq cram bams filter_mutect2
