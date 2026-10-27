#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -l mem=8G
#$ -pe smp 2
# One task per patient. Replaces pyclone_counts.sh, pileup.sh,
# build_pyclone_input.sh, pyclone_run.sh and summarize_pyclone.sh, which had
# the patient/sample names typed in (and mismatched: 037 vs 063, BL vs BSL)
# and called a summary script name that did not exist.
# usage: pyclone.job.sh <list of patients>
source "$PIPELINE_DIR/lib/common.sh"
p=$(task_line < "$1")
dir=$PYCLONE_DIR/$p
man=$dir/manifest.tsv

# conda env first (bcftools, alleleCounter, pyclone-vi on PATH), then R from
# the module system. TODO: if R and the env's libraries clash, install
# r-data.table + r-ggplot2 into the env instead.
use_conda pyclone-vi
module unload gcc-libs 2> /dev/null || true
module unload mpi 2> /dev/null || true
module load "$R_MODULE"

# 1. union of PASS SNVs across this patient's tumours -> <p>.union_sites.tsv
Rscript "$PIPELINE_DIR/bin/build_pyclone_input.R" --emit-sites "$man" "$dir"
sites=$dir/$p.union_sites.tsv
cut -f1,2 "$sites" > "$dir/$p.loci.txt"

# 2. forced allele counts at every union site in every tumour sample
mkdir -p "$dir/counts"
tail -n +2 "$man" | while IFS=$'\t' read -r _ s _ cram counts _ _ _; do
    alleleCounter -l "$dir/$p.loci.txt" -r "$REF_FASTA" -m 20 -q 20 \
        -b "$cram" -o "$dir/counts/$s.allelecount.txt"
    # alleleCounter: #CHR POS Count_A Count_C Count_G Count_T Good_depth
    awk 'BEGIN { OFS = "\t"; b["A"] = 3; b["C"] = 4; b["G"] = 5; b["T"] = 6 }
         NR == FNR { ref[$1 ":" $2] = $3; alt[$1 ":" $2] = $4; next }
         /^#/ { next }
         { k = $1 ":" $2; if (k in ref) print $1, $2, ref[k], alt[k], $(b[ref[k]]), $(b[alt[k]]) }' \
        "$sites" "$dir/counts/$s.allelecount.txt" > "$counts"
    log "$s: $(wc -l < "$counts") sites counted"
done

# 3. PyClone-VI input table -> <p>.pyclone_vi.tsv
Rscript "$PIPELINE_DIR/bin/build_pyclone_input.R" "$man" "$dir" vi

# 4. fit
pyclone-vi fit -i "$dir/$p.pyclone_vi.tsv" -o "$dir/$p.h5" \
    -c 40 -d beta-binomial -r 10 --num-threads 2
pyclone-vi write-results-file -i "$dir/$p.h5" -o "$dir/$p.results.tsv"

# 5. cluster table + timepoint CCF plot
Rscript "$PIPELINE_DIR/bin/summarize_pyclone.R" "$dir/$p.results.tsv" "$dir"
