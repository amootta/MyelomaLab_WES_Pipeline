#!/bin/bash -l
#$ -l h_rt=2:00:00
#$ -l mem=8G
#$ -pe smp 1
# mgp1000 union consensus TSV -> filter_consensus_calls.py -> consensus_tsv_to_vcf.py,
# then file each VCF under pvacseq_input/<patient>/ by its tumour sample name.
# Replaces run_consensus_conversion.sh.
#
# mgp1000's output folder is shared by every batch, so only pairs whose tumour
# is in this batch's samples.tsv are taken (see sheet_tumour in lib/common.sh).
source "$PIPELINE_DIR/lib/common.sh"
module load "$PYTHON_MODULE"

tables=$CONSENSUS_DIR/mgp1000_tables
filtered=$CONSENSUS_DIR/filtered
vcfs=$CONSENSUS_DIR/vcf
by_patient=$CONSENSUS_DIR/pvacseq_input
mkdir -p "$tables" "$filtered" "$vcfs" "$by_patient"
# rebuilt from scratch every run, so a changed sample sheet leaves nothing stale
find "$tables" "$by_patient" -type l -delete
find "$filtered" "$vcfs" -type f -delete

kept=0 skipped=0
for f in $MGP1000_CONSENSUS_GLOB; do
    [[ -e $f ]] || continue
    pair=$(basename "$f")
    pair=${pair%%.hq.union.*}
    if [[ -n $(pair_tumour "$pair") ]]; then
        ln -s "$f" "$tables/$(basename "$f")"
        kept=$((kept + 1))
    else
        skipped=$((skipped + 1))
    fi
done
log "$kept pair(s) from this batch, $skipped other pair(s) in $MGP1000_SOMATIC_DIR/consensus skipped"
((kept > 0)) || die "no mgp1000 consensus tables match the tumours in $SAMPLES"

python3 "$PIPELINE_DIR/bin/filter_consensus_calls.py" \
    --input "$tables/*.hq.union.consensus.somatic.snv.txt.gz" \
    --min-callers "$CONSENSUS_MIN_CALLERS" \
    --outdir "$filtered"

python3 "$PIPELINE_DIR/bin/consensus_tsv_to_vcf.py" \
    --input "$filtered/*.callers_ge${CONSENSUS_MIN_CALLERS}.txt*" \
    --outdir "$vcfs" \
    --reference-fai "$REF_FAI"

declare -A seen
for v in "$vcfs"/*.vcf; do
    name=$(grep -m1 '^#CHROM' "$v" | cut -f10)
    tumour=$(sheet_tumour "$name")
    [[ -n $tumour ]] || { warn "$v: tumour '$name' is not in $SAMPLES, not filed"; continue; }
    if [[ -n ${seen[$tumour]:-} ]]; then
        warn "$tumour has two consensus VCFs (${seen[$tumour]} and $(basename "$v")) -- an old and a new mgp1000 run? Both filed; remove one from mgp1000's output."
    fi
    seen[$tumour]=$(basename "$v")
    patient=$(patient_of "$tumour")
    mkdir -p "$by_patient/$patient"
    ln -sf "$v" "$by_patient/$patient/$(basename "$v")"
    log "$name -> $tumour ($patient)"
done
