#!/bin/bash -l
#$ -l h_rt=1:00:00
#$ -l mem=4G
#$ -pe smp 1
# T1K KIR calls vs ground-truth typing -> workbook.
source "$PIPELINE_DIR/lib/common.sh"
module load "$PYTHON_MODULE"

python3 "$PIPELINE_DIR/bin/kir_concordance.py" \
    --t1k-dir "$KIR_DIR/T1K" \
    --filename-pattern 'cos_(\d+)_(\d+)_KIR_genotype' \
    ${KIR_TRUTH:+--truth-csv "$KIR_TRUTH"} \
    --out "$KIR_DIR/KIR_concordance.xlsx"
# TODO: KIRCLE output is not in the concordance report yet (kir_concordance.py
#       only reads T1K).
