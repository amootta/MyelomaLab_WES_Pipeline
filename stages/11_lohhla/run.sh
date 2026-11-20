#!/bin/bash
# Stage 11 -- LOHHLA: HLA class I loss of heterozygosity, per patient (all of a
# patient's tumours against its normal in one run). Details: docs/LOHHLA.md.
#
#   ./wes.sh lohhla setup             once per REF_DIR: clone LOHHLA, patch it for hg38 +
#                                     picard.jar, fetch the IMGT HLA fasta/hla.dat
#   ./wes.sh lohhla inputs            per patient: bam links, HLA alleles, purity/ploidy
#                                                                   -> LOHHLA/<P>/input/
#   ./wes.sh lohhla run [patient..]   one SGE task per patient (all, or those named)
#                                                                   -> LOHHLA/<P>/
#   ./wes.sh lohhla summary [opts]    LOH call per tumour x gene, all finished patients
#                                                                   -> LOHHLA/LOHHLA_summary{,_wide}.tsv
#
# Needs: stage 01 bams + ASCAT, stage 04 `hla collect`, and the `lohhla` conda
# env (conda env create -f envs/lohhla.yml).
source "$(dirname "$0")/../../lib/common.sh"
here=$PIPELINE_DIR/stages/11_lohhla
source "$here/lohhla_env.sh"

LOHHLA_GIT=https://bitbucket.org/mcgranahanlab/lohhla.git
IMGT_URL=https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla

step=${1:-}
shift || true
ARGS=("$@")   # extra arguments: patients for `run`, options for `summary`

step_setup() {
    use_conda "$LOHHLA_CONDA_ENV"
    mkdir -p "$LOHHLA_HOME/imgt"
    if [[ ! -d $LOHHLA_HOME/src/.git ]]; then
        log "cloning $LOHHLA_GIT"
        git clone "$LOHHLA_GIT" "$LOHHLA_HOME/src"
    fi
    log "LOHHLA commit: $(git -C "$LOHHLA_HOME/src" rev-parse --short HEAD)"

    # 1. The script's three HLA-A/B/C read windows are hg19 coordinates. Moved
    #    to GRCh38 (+32,223 bp across this part of chr6). Reads are also
    #    "fished" from the whole bam by HLA k-mers, so these windows only add to that.
    # 2. It calls separate Picard jars (<gatkDir>/SamToFastq.jar ...); current
    #    Picard is one picard.jar, so call the tool by name instead.
    sed -e 's/chr6:29909037-29913661/chr6:29941260-29945884/' \
        -e 's/chr6:31321649-31324964/chr6:31353872-31357187/' \
        -e 's/chr6:31236526-31239869/chr6:31268749-31272092/' \
        -e 's#/SamToFastq\.jar#/picard.jar SamToFastq#g' \
        -e 's#/SortSam\.jar#/picard.jar SortSam#g' \
        -e 's#/FilterSamReads\.jar#/picard.jar FilterSamReads#g' \
        "$LOHHLA_HOME/src/LOHHLAscript.R" > "$LOHHLA_SCRIPT"
    local s
    for s in chr6:29941260 chr6:31353872 chr6:31268749 "picard.jar SamToFastq" "picard.jar SortSam" "picard.jar FilterSamReads"; do
        grep -qF "$s" "$LOHHLA_SCRIPT" || die "patch did not apply ($s) -- has LOHHLAscript.R changed? Check $LOHHLA_SCRIPT by hand"
    done
    log "wrote $LOHHLA_SCRIPT (hg38 windows, picard.jar)"

    # IMGT/HLA genomic sequences, class I A/B/C only, renamed the way LOHHLA
    # expects: ">HLA:HLA00001 A*01:01:01:01 3503 bp" -> ">hla_a_01_01_01_01"
    local imgt=$LOHHLA_HOME/imgt
    if [[ ! -s $imgt/hla_gen.fasta ]]; then
        log "downloading IMGT hla_gen.fasta and hla.dat"
        curl -fsSL -o "$imgt/hla_gen.fasta" "$IMGT_URL/hla_gen.fasta"
        curl -fsSL -o "$imgt/hla.dat" "$IMGT_URL/hla.dat"
        curl -fsSL "$IMGT_URL/Allelelist.txt" 2> /dev/null | grep -m1 -i version > "$imgt/VERSION" || true
        echo "downloaded $(date +%F)" >> "$imgt/VERSION"
    fi
    if [[ $LOHHLA_HLA_FASTA == "$imgt/hla_abc_gen.fasta" && ! -s $LOHHLA_HLA_FASTA ]]; then
        awk '/^>/ { split($2, a, "*"); keep = (a[1] == "A" || a[1] == "B" || a[1] == "C")
                    if (keep) { n = tolower("hla_" a[1] "_" a[2]); gsub(":", "_", n); print ">" n }
                    next }
             keep' "$imgt/hla_gen.fasta" > "$LOHHLA_HLA_FASTA"
        samtools faidx "$LOHHLA_HLA_FASTA"
    fi
    log "HLA fasta: $LOHHLA_HLA_FASTA ($(grep -c '>' "$LOHHLA_HLA_FASTA") alleles); $(head -1 "$imgt/VERSION" 2> /dev/null)"
    novo_setup
    log "novoalign: $novo_dir/novoalign ($("$novo_dir/novoalign" 2>&1 | grep -m1 -i 'novoalign' || true))"
    ls "$CONDA_PREFIX"/share/picard-*/picard.jar > /dev/null 2>&1 || warn "picard.jar not found under $CONDA_PREFIX/share"
}

step_inputs() {
    need "$LOHHLA_HLA_FASTA" "$HLA_DIR/calls/$LOHHLA_HLA_SOURCE"
    local manifest=$SAMPLESHEET_DIR/lohhla_manifest.tsv p t n
    mkdir -p "$SAMPLESHEET_DIR" "$LOHHLA_DIR"
    : > "$manifest"
    declare -A seen=()
    while IFS=$'\t' read -r p t n; do
        if [[ -z ${seen[$p]:-} ]]; then
            need "$(recal_bam "$n")"
            printf '%s\t%s\tnormal\t%s\t\n' "$p" "$n" "$(recal_bam "$n")" >> "$manifest"
            seen[$p]=1
        fi
        need "$(recal_bam "$t")"
        [[ -e $(ascat_pp "$t" "$n") ]] || warn "$t: no ASCAT purity/ploidy -- LOHHLA reports allelic imbalance only"
        printf '%s\t%s\ttumour\t%s\t%s\n' "$p" "$t" "$(recal_bam "$t")" "$(ascat_pp "$t" "$n")" >> "$manifest"
    done < <(pairs)
    log "wrote $manifest"
    mod_load "$PYTHON_MODULE"
    python3 -c 'import pandas' 2> /dev/null || die "python3 ($(command -v python3)) has no pandas -- check PYTHON_MODULE=$PYTHON_MODULE"
    python3 "$PIPELINE_DIR/bin/build_lohhla_inputs.py" \
        --manifest "$manifest" \
        --calls-dir "$HLA_DIR/calls" \
        --source "$LOHHLA_HLA_SOURCE" \
        --hla-fasta "$LOHHLA_HLA_FASTA" \
        --outdir "$LOHHLA_DIR"
    echo "Check $LOHHLA_DIR/<patient>/input/{hla_alleles,solutions}.txt before 'wes.sh lohhla run'." >&2
}

step_run() {
    local list=$SAMPLESHEET_DIR/lohhla_patients.txt p
    need "$LOHHLA_SCRIPT" "$LOHHLA_HLA_FASTA"
    mkdir -p "$SAMPLESHEET_DIR"
    : > "$list"
    local want=("${ARGS[@]}")
    ((${#want[@]})) || mapfile -t want < <(patients)
    for p in "${want[@]}"; do
        if [[ -s $LOHHLA_DIR/$p/input/hla_alleles.txt ]]; then
            echo "$p" >> "$list"
        else
            warn "$p: no inputs (run 'wes.sh lohhla inputs') -- skipped"
        fi
    done
    submit lohhla "$(wc -l < "$list")" "$here/lohhla.job.sh" "$list"
}

# Extra arguments go to summarize_lohhla.py (e.g. --min-sites 10 --pval 0.05).
step_summary() {
    compgen -G "$LOHHLA_DIR/*/*HLAlossPrediction_CI*" > /dev/null || die "no finished LOHHLA runs in $LOHHLA_DIR yet"
    mod_load "$PYTHON_MODULE"
    python3 "$PIPELINE_DIR/bin/summarize_lohhla.py" --lohhla-dir "$LOHHLA_DIR" --out "$LOHHLA_DIR/LOHHLA_summary" "${ARGS[@]}"
}

run_steps "usage: wes.sh lohhla <setup|inputs|run [patient...]|summary [options]>" "$step" setup inputs run summary
