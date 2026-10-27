#!/bin/bash -l
#$ -l h_rt=40:00:00
#$ -l mem=16G
#$ -pe smp 1
# Nextflow head job for nf-core/sarek. Submitted by stages/01_sarek/run.sh.
# usage: sarek.job.sh <fastq|cram> <samplesheet.csv>
source "$PIPELINE_DIR/lib/common.sh"
mode=$1 sheet=$2
nextflow_env

mode_args=(--trim_fastq --save_mapped)      # fastq: full run from mapping
[[ $mode == cram ]] && mode_args=(--step variant_calling)

mkdir -p "$SAREK_DIR"
cd "$SAREK_DIR"
nextflow run nf-core/sarek ${SAREK_VERSION:+-r "$SAREK_VERSION"} \
    -profile ucl_myriad \
    -c "$PIPELINE_DIR/config/nextflow/myriad_extra.config" \
    -resume \
    --input "$sheet" \
    --outdir "$SAREK_DIR" \
    --wes \
    "${mode_args[@]}" \
    --tools "$SAREK_TOOLS" \
    --igenomes_ignore \
    --fasta "$REF_FASTA" \
    --fasta_fai "$REF_FAI" \
    --dbsnp "$DBSNP" \
    --known_indels "$KNOWN_INDELS" \
    --ascat_genome hg38 \
    --ascat_loci "$ASCAT_LOCI" \
    --ascat_alleles "$ASCAT_ALLELES" \
    --ascat_loci_gc "$ASCAT_LOCI_GC" \
    --ascat_loci_rt "$ASCAT_LOCI_RT" \
    --mappability "$MAPPABILITY" \
    --chr_dir "$CONTROLFREEC_CHR_DIR" \
    --email "$EMAIL" \
    --max_cpus 8
