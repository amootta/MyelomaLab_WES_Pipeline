#!/bin/bash -l
#$ -l h_rt=2:00:00
#$ -l mem=8G
#$ -pe smp 1
# mgp1000 union consensus TSV -> filter_consensus_calls.py -> consensus_tsv_to_vcf.py,
# then file each VCF under pvacseq_input/<patient>/ by its tumour sample name.
# Replaces run_consensus_conversion.sh.
source "$PIPELINE_DIR/lib/common.sh"
module load "$PYTHON_MODULE"

filtered=$CONSENSUS_DIR/filtered
vcfs=$CONSENSUS_DIR/vcf
by_patient=$CONSENSUS_DIR/pvacseq_input
mkdir -p "$filtered" "$vcfs" "$by_patient"

python3 "$PIPELINE_DIR/bin/filter_consensus_calls.py" \
    --input "$MGP1000_CONSENSUS_GLOB" \
    --min-callers "$CONSENSUS_MIN_CALLERS" \
    --outdir "$filtered"

python3 "$PIPELINE_DIR/bin/consensus_tsv_to_vcf.py" \
    --input "$filtered/*.callers_ge${CONSENSUS_MIN_CALLERS}.txt*" \
    --outdir "$vcfs" \
    --reference-fai "$REF_FAI"

# Tumour names in mgp1000 output can carry lane/library suffixes
# (COS_001_063_BL_EKDN...L3), so match on the longest patient ID prefix.
mapfile -t pts < <(patients | awk '{ print length, $0 }' | sort -rn | cut -d' ' -f2-)
for v in "$vcfs"/*.vcf; do
    tumour=$(grep -m1 '^#CHROM' "$v" | cut -f10)
    patient=""
    for p in "${pts[@]}"; do
        [[ $tumour == "$p"* ]] && { patient=$p; break; }
    done
    if [[ -z $patient ]]; then
        warn "$v: tumour '$tumour' matches no patient in $SAMPLES, not filed"
        continue
    fi
    mkdir -p "$by_patient/$patient"
    ln -sf "$v" "$by_patient/$patient/$(basename "$v")"
    log "$tumour -> $patient"
done
