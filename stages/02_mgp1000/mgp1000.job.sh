#!/bin/bash -l
#$ -l h_rt=48:00:00
#$ -l mem=16G
#$ -pe smp 2
# Nextflow head job for one mgp1000 workflow (it submits its own SGE jobs via
# --executor sge). From mgp1000_WES.sh. The head job only orchestrates, but it
# has to outlive the whole run; if it hits the wall-time limit, rerun the same
# wes.sh command and -resume carries on.
# usage: mgp1000.job.sh <somatic|germline> <sample sheet>
source "$PIPELINE_DIR/lib/common.sh"
workflow=$1 sheet=$2

for m in $MGP1000_APPTAINER_MODULES; do module load "$m"; done
module load samtools
module load "$JAVA_MODULE"
export NXF_VER=$MGP1000_NXF_VER
export NXF_SINGULARITY_CACHEDIR=$MGP1000_SINGULARITY_CACHE

case $workflow in
    somatic)  tools=$MGP1000_SOMATIC_TOOLS ;;
    germline) tools=$MGP1000_GERMLINE_TOOLS ;;
    *) die "unknown mgp1000 workflow '$workflow'" ;;
esac

# The install in MGP1000_DIR is shared by every batch; launching from this
# batch's MGP1000_RUN_DIR keeps work/, .nextflow/ and output/ per batch.
mkdir -p "$MGP1000_RUN_DIR"
cd "$MGP1000_RUN_DIR"
# shellcheck disable=SC2086  # $tools is a list of flags
nextflow run "$MGP1000_DIR/$workflow.nf" \
    --run_id "${MGP1000_RUN_ID}_$workflow" \
    --seq_protocol WES \
    --sample_sheet "$sheet" \
    --input_dir "$MGP1000_INPUT_DIR" \
    --output_dir "$MGP1000_RUN_DIR/output" \
    --email "$EMAIL" \
    --executor sge \
    -profile "$workflow" \
    $tools \
    -resume

# --- not wired into wes.sh ---------------------------------------------------
# Preprocessing from FASTQs (instead of sarek's bams), once per sample, before
# germline/somatic. FASTQs must be named *_R1/_R2 and end in .fq.gz:
#   nextflow run preprocessing.nf --run_id <sample> --input_format fastq \
#     --input_dir <fastq dir> --email "$EMAIL" -profile preprocessing \
#     --executor sge --seq_protocol WES --cpus 8 --output_dir .
