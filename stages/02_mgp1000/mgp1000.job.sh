#!/bin/bash -l
#$ -l h_rt=48:00:00
#$ -l mem=16G
#$ -pe smp 1
# Nextflow head job for mgp1000's somatic workflow.
# usage: mgp1000.job.sh <mgp1000.csv>
source "$PIPELINE_DIR/lib/common.sh"
sheet=$1
nextflow_env
export NXF_VER=$MGP1000_NXF_VER   # mgp1000 was written for nextflow 22.10.x

cd "$MGP1000_DIR"
# TODO: replace this block with the command from your mgp1000_WES.sh (never
#       committed to the repo). Keep --sample_sheet "$sheet" and the cohort
#       name so stage 03 finds output/somatic/.
die "mgp1000 command not filled in yet -- edit stages/02_mgp1000/mgp1000.job.sh"
# nextflow run somatic.nf -resume \
#     --sample_sheet "$sheet" \
#     --cohort_name "$BATCH_NAME" \
#     --email "$EMAIL" \
#     -profile somatic
