#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -l mem=16G
#$ -pe smp 2
# Nextflow head job for mgp1000's somatic workflow (it submits its own SGE
# jobs via --executor sge). From mgp1000_WES.sh.
# usage: mgp1000.job.sh <mgp1000_somatic.csv>
source "$PIPELINE_DIR/lib/common.sh"
sheet=$1

for m in $MGP1000_APPTAINER_MODULES; do module load "$m"; done
module load samtools
module load "$JAVA_MODULE"
export NXF_VER=$MGP1000_NXF_VER
export NXF_SINGULARITY_CACHEDIR=$MGP1000_SINGULARITY_CACHE

cd "$MGP1000_DIR"
# shellcheck disable=SC2086  # MGP1000_SOMATIC_TOOLS is a list of flags
nextflow run ./somatic.nf \
    --run_id "$MGP1000_RUN_ID" \
    --seq_protocol WES \
    --sample_sheet "$sheet" \
    --input_dir "$MGP1000_INPUT_DIR" \
    --email "$EMAIL" \
    --executor sge \
    -profile somatic \
    $MGP1000_SOMATIC_TOOLS \
    -resume

# --- alternatives from mgp1000_WES.sh, not wired into wes.sh ----------------
# Germline, per sample:
#   nextflow run ./germline.nf --run_id <sample> --seq_protocol WES \
#     --sample_sheet <sheet> --input_dir <dir with bams> --email "$EMAIL" \
#     -profile germline --executor sge
#
# Preprocessing from FASTQs (instead of sarek's bams), once per sample, before
# germline/somatic. FASTQs must be named *_R1/_R2 and end in .fq.gz:
#   nextflow run preprocessing.nf --run_id <sample> --input_format fastq \
#     --input_dir <fastq dir> --email "$EMAIL" -profile preprocessing \
#     --executor sge --seq_protocol WES --cpus 8 --output_dir .
