#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -l mem=6G
#$ -pe smp 5
# One task per patient: LOHHLA on every tumour bam in <P>/input/bams against the normal.
# usage: lohhla.job.sh <patients.txt>
# Test one patient interactively, without SGE:
#   TASK=1 PIPELINE_DIR=$PWD bash stages/11_lohhla/lohhla.job.sh $BATCH_DIR/samplesheets/lohhla_patients.txt
source "$PIPELINE_DIR/lib/common.sh"
source "$PIPELINE_DIR/stages/11_lohhla/lohhla_env.sh"
p=$(task_line < "$1")
use_conda "$LOHHLA_CONDA_ENV"

in=$LOHHLA_DIR/$p/input
out=$LOHHLA_DIR/$p
normal=$(< "$in/normal.txt")
picard_dir=$(dirname "$(ls "$CONDA_PREFIX"/share/picard-*/picard.jar | head -1)")
novo_setup
need "$in/bams/$normal.bam" "$in/hla_alleles.txt" "$LOHHLA_SCRIPT" "$LOHHLA_HLA_FASTA" "$picard_dir/picard.jar"

# LOHHLA runs every tool through R's system() and ignores failures: a missing
# jellyfish just leaves an empty k-mer list, the whole bam is grepped for
# nothing, and it dies later in combine.fastqs. So check the tools first.
for tool in jellyfish samtools bedtools java novoindex; do
    command -v "$tool" > /dev/null || die "$tool not on PATH"
done
jellyfish --version 2>&1 | grep -q 'jellyfish 2' ||
    die "LOHHLA needs jellyfish 2 (writes mer_counts.jf); got: $(jellyfish --version 2>&1 | head -1)"
log "tools: $(jellyfish --version 2>&1 | head -1), $(samtools --version | head -1), bedtools $(bedtools --version | cut -d' ' -f2), novoalign in $novo_dir"

# Integer copy number only if at least one tumour has an ASCAT solution.
cn=FALSE
(($(wc -l < "$in/solutions.txt") > 1)) && cn=$in/solutions.txt

log "LOHHLA $p: normal $normal, tumours: $(ls "$in/bams" | grep '\.bam$' | grep -vx "$normal.bam" | tr '\n' ' ')"
cd "$out"
# BAMDir and normalBAMfile must agree as strings: LOHHLA spots the normal by
# comparing "<BAMDir>/<file>" with --normalBAMfile.
Rscript "$LOHHLA_SCRIPT" \
    --patientId "$p" \
    --outputDir "$out" \
    --normalBAMfile "$in/bams/$normal.bam" \
    --BAMDir "$in/bams" \
    --hlaPath "$in/hla_alleles.txt" \
    --HLAfastaLoc "$LOHHLA_HLA_FASTA" \
    --HLAexonLoc "$LOHHLA_HLA_DAT" \
    --CopyNumLoc "$cn" \
    --minCoverageFilter "$LOHHLA_MIN_COV" \
    --mappingStep TRUE \
    --fishingStep TRUE \
    --plottingStep "$LOHHLA_PLOTS" \
    --coverageStep TRUE \
    --cleanUp "$LOHHLA_CLEANUP" \
    --gatkDir "$picard_dir" \
    --novoDir "$novo_dir"

ls "$out"/*HLAlossPrediction* > /dev/null 2>&1 || die "no HLAlossPrediction table in $out -- see the running.hla.loh.exome@*_log.txt there"
log "done: $(ls "$out"/*HLAlossPrediction*)"
