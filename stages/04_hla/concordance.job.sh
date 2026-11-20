#!/bin/bash -l
#$ -l h_rt=1:00:00
#$ -l mem=4G
#$ -pe smp 1
# OptiType / HLA-HD / T1K vs ground-truth typing -> colour-coded workbook.
source "$PIPELINE_DIR/lib/common.sh"
module load "$PYTHON_MODULE"

calls=$HLA_DIR/calls
truth_flag=--truth-csv
[[ $HLA_TRUTH == *.xlsx || $HLA_TRUTH == *.xls ]] && truth_flag=--truth-xlsx
python3 "$PIPELINE_DIR/bin/hla_concordance.py" \
    --optitype-dir "$calls/optitype" \
    --hlahd-dir "$calls/hlahd" \
    --t1k-dir "$calls/t1k" \
    "$truth_flag" "$HLA_TRUTH" \
    --out "$HLA_DIR/HLA_typing_comparison.xlsx"
