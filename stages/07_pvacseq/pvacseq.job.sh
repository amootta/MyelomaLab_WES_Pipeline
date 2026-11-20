#!/bin/bash -l
#$ -l h_rt=40:00:00
#$ -l mem=16G
#$ -pe smp 1
# One task per patient: nextflow-pvacseq head job on that patient's tumour VCFs.
# usage: pvacseq.job.sh <list of patients>
source "$PIPELINE_DIR/lib/common.sh"
p=$(task_line < "$1")
nextflow_env   # sets TMPDIR on Scratch -- needed for the pvactools image pull

out=$PVACSEQ_DIR/$p
mkdir -p "$out"
cd "$out"
# A cache folder already in place saves the ~20 GB VEP download in every
# patient's run; without one the pipeline downloads it itself.
vep_cache=()
[[ -n ${PVACSEQ_VEP_CACHE:-} && -d $PVACSEQ_VEP_CACHE ]] && vep_cache=(--vep_cache "$PVACSEQ_VEP_CACHE" --vep_cache_version "$PVACSEQ_VEP_CACHE_VERSION")
nextflow run "$PVACSEQ_PIPELINE" \
    -profile "$PVACSEQ_PROFILE" \
    -c "$PIPELINE_DIR/config/nextflow/myriad_pvacseq.config" \
    -resume \
    --max_time 40h \
    --input "$CONSENSUS_DIR/pvacseq_input/$p/" \
    --hla_csv "$PVACSEQ_DIR/hla/$p.csv" \
    --fasta "$REF_FASTA" \
    --outdir "$out" \
    --pvacseq_iedb "$PVACSEQ_IEDB" \
    --vep_plugins "$PVACSEQ_VEP_PLUGINS" \
    "${vep_cache[@]}" \
    --pvacseq_algorithm "$PVACSEQ_ALGORITHMS"
