#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -l mem=8G
#$ -pe smp 1
# One task per tumour/normal pair.
# Uses the pre-BQSR mapped bams (sarek --save_mapped), as the old
# igcaller_wes.sh did; falls back to the recal bams when those are absent.
# usage: igcaller.job.sh <pairs.tsv: patient tumour normal>
source "$PIPELINE_DIR/lib/common.sh"
IFS=$'\t' read -r p t n < <(task_line < "$1")
module load python3

out=$IGCALLER_DIR/$t
mkdir -p "$out"
cd "$IGCALLER_HOME"
python3 IgCaller \
    -V hg38 \
    -I IgCaller_reference_files/ \
    -C ucsc \
    -T "$(mapped_bam "$t")" \
    -N "$(mapped_bam "$n")" \
    -R "$REF_FASTA" \
    -o "$out/" \
    -seq wes \
    -seqDepth high \
    -pm paired \
    -ptsam "$IGCALLER_SAMTOOLS"
