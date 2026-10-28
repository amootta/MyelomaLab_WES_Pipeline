#!/bin/bash
# Stage 02 -- MGP1000 on the recal bams from stage 01. mgp1000 only reads
# bams, so run `wes.sh sarek bams` first.
#
#   ./wes.sh mgp1000 somatic    per pair: Mutect2 + Strelka2 + VarScan consensus,
#                               Manta, Delly, IgCaller   -> output/somatic/
#   ./wes.sh mgp1000 germline   per normal: DeepVariant germline SNVs/indels and
#                               fastNGSadmix ancestry    -> output/germline/
#
# Starting from FASTQs instead (mgp1000's own preprocessing.nf, not sarek's
# bams): see the commented preprocessing.nf command in mgp1000.job.sh.
source "$(dirname "$0")/../../lib/common.sh"
here=$PIPELINE_DIR/stages/02_mgp1000

# Put each bam + index into one --input_dir as <sample>.final.bam, mgp1000's
# own naming: somatic.nf takes the ID as everything before the first ".", and
# germline.nf strips ".final...bam", so both name outputs after the sample.
# Hard links, not symlinks: mgp1000's containers only see the bound input
# folder, so a symlink pointing elsewhere on Scratch would dangle inside them.
# Falls back to a copy across filesystems.
stage_bam() {  # stage_bam <sample>
    local src dst idx
    src=$(recal_bam "$1") dst=$MGP1000_INPUT_DIR/$1.final.bam
    need "$src"
    [[ $src -ef $dst ]] || ln -f "$src" "$dst" 2> /dev/null || cp "$src" "$dst"
    for idx in "$src.bai" "${src%.bam}.bai"; do
        if [[ -e $idx ]]; then
            [[ $idx -ef $dst.bai ]] || ln -f "$idx" "$dst.bai" 2> /dev/null || cp "$idx" "$dst.bai"
            return
        fi
    done
    die "no index for $src"
}

step_somatic() {
    local sheet=$SAMPLESHEET_DIR/mgp1000_somatic.csv t n
    mkdir -p "$SAMPLESHEET_DIR" "$MGP1000_INPUT_DIR"
    echo "normal,tumor" > "$sheet"
    while IFS=$'\t' read -r _ t n; do
        stage_bam "$t"
        stage_bam "$n"
        echo "$n.final.bam,$t.final.bam" >> "$sheet"
    done < <(pairs)
    log "wrote $sheet; bams in $MGP1000_INPUT_DIR"
    submit mgp1000_somatic 1 "$here/mgp1000.job.sh" somatic "$sheet"
}

# germline.nf reads only the normal column, so one row per normal (a normal
# listed once per tumour would be called several times).
step_germline() {
    local sheet=$SAMPLESHEET_DIR/mgp1000_germline.csv p n t f
    mkdir -p "$SAMPLESHEET_DIR" "$MGP1000_INPUT_DIR"
    echo "normal,tumor" > "$sheet"
    while IFS=$'\t' read -r p n; do
        stage_bam "$n"
        t=$(tumours | awk -F'\t' -v p="$p" '$1 == p { print $2; exit }')
        echo "$n.final.bam,${t:+$t.final.bam}" >> "$sheet"
    done < <(normals)
    for f in fastNGSadmix_markers_refPanel.hg38.txt fastNGSadmix_markers_nInd.hg38.txt fastNGSadmix_markers.hg38.sites; do
        [[ -e $MGP1000_DIR/references/hg38/$f ]] ||
            warn "ancestry reference $MGP1000_DIR/references/hg38/$f not found -- fastNGSadmix will fail"
    done
    log "wrote $sheet ($(($(wc -l < "$sheet") - 1)) normals); bams in $MGP1000_INPUT_DIR"
    submit mgp1000_germline 1 "$here/mgp1000.job.sh" germline "$sheet"
}

run_steps "usage: wes.sh mgp1000 <somatic|germline>" "${1:-}" somatic germline
