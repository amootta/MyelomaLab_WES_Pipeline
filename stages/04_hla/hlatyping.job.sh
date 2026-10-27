#!/bin/bash -l
#$ -l h_rt=40:00:00
#$ -l mem=16G
#$ -pe smp 1
# Nextflow head job for the patched nf-core/hlatyping checkout (OptiType + HLA-HD).
# From nextflow_HLA.sh. Runs the local checkout in HLATYPING_DIR, not the
# released pipeline.
# usage: hlatyping.job.sh <hlatyping.csv>
source "$PIPELINE_DIR/lib/common.sh"
sheet=$1
nextflow_env

mkdir -p "$HLA_DIR/hlatyping"
cd "$HLATYPING_DIR"
nextflow run . \
    -profile ucl_myriad \
    -c "$PIPELINE_DIR/config/nextflow/hlatyping_fix.config" \
    -resume \
    --input "$sheet" \
    --outdir "$HLA_DIR/hlatyping" \
    --tools optitype,hlahd \
    --hlahd_path "$HLAHD_TARBALL" \
    --hlahd_update_reference_dict true
