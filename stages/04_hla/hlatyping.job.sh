#!/bin/bash -l
#$ -l h_rt=40:00:00
#$ -l mem=16G
#$ -pe smp 1
# Nextflow head job for nf-core/hlatyping (OptiType + HLA-HD).
# usage: hlatyping.job.sh <hlatyping.csv>
source "$PIPELINE_DIR/lib/common.sh"
sheet=$1
nextflow_env

mkdir -p "$HLA_DIR/hlatyping"
cd "$HLA_DIR/hlatyping"
# TODO: nextflow_HLA.sh was never committed. Check the release and how HLA-HD
#       was switched on (the hlahd.1.7.1.tar.gz from the RDSS), and pin
#       HLATYPING_VERSION in the config.
nextflow run nf-core/hlatyping ${HLATYPING_VERSION:+-r "$HLATYPING_VERSION"} \
    -profile ucl_myriad \
    -c "$PIPELINE_DIR/config/nextflow/hlatyping_fix.config" \
    -resume \
    --input "$sheet" \
    --outdir "$HLA_DIR/hlatyping" \
    --genome GRCh38 \
    --email "$EMAIL"
