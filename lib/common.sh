# =============================================================================
# common.sh -- sourced by every launcher (stages/*/run.sh) and every SGE job
# (stages/*/*.job.sh). Loads the config, reads config/samples.tsv, and knows
# where each stage puts its files, so no script has to hardcode a path.
#
# Launchers run on the login node, check that their inputs exist, then qsub
# array jobs. SGE copies job scripts to a spool dir, so jobs can't find the
# repo from their own location -- `submit` passes PIPELINE_DIR /
# PIPELINE_CONFIG / SAMPLES through with qsub -v.
# =============================================================================

if [[ -z ${PIPELINE_DIR:-} ]]; then
    PIPELINE_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fi
PIPELINE_CONFIG=${PIPELINE_CONFIG:-$PIPELINE_DIR/config/pipeline.config}
_samples_override=${SAMPLES:-}
[[ -f $PIPELINE_CONFIG ]] || { echo "ERROR: config not found: $PIPELINE_CONFIG" >&2; exit 1; }
# shellcheck source=../config/pipeline.config
source "$PIPELINE_CONFIG"
[[ -n $_samples_override ]] && SAMPLES=$_samples_override
[[ -f $SAMPLES ]] || { echo "ERROR: sample sheet not found: $SAMPLES" >&2; exit 1; }
export PIPELINE_DIR PIPELINE_CONFIG SAMPLES

# Not -u: Myriad's `module` and conda's activate scripts reference unset vars.
set -eo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] $*" >&2; }
warn() { echo "WARNING: $*" >&2; }

# need <file>... : stop with a clear message if any input is missing.
need() {
    local f missing=0
    for f in "$@"; do
        [[ -e $f ]] || { echo "  missing: $f" >&2; missing=1; }
    done
    ((missing == 0)) || die "required input(s) missing -- has the previous stage finished?"
}

# have <file>... : like need, but only reports; returns 1 if anything is missing.
have() {
    local f
    for f in "$@"; do [[ -e $f ]] || return 1; done
}

# -----------------------------------------------------------------------------
# samples.tsv  (tab-separated, header row; one row per sequencing LANE)
#   patient sample status(0=normal,1=tumour) sex lane fastq_1 fastq_2 cram bam
# cram/bam are optional overrides for samples that did not come out of this
# batch's sarek run (merged crams, germline bams from elsewhere, ...).
# -----------------------------------------------------------------------------

# sample_rows: one line per SAMPLE (lanes collapsed): patient sample status sex
sample_rows() {
    awk -F'\t' '
        NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
        /^#/ || NF == 0 { next }
        !seen[$c["sample"]]++ { print $c["patient"] "\t" $c["sample"] "\t" $c["status"] "\t" $c["sex"] }
    ' "$SAMPLES"
}

# sample_field <sample> <column>: first non-empty value of a column for a sample
sample_field() {
    awk -F'\t' -v s="$1" -v col="$2" '
        NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
        $c["sample"] == s && $c[col] != "" { print $c[col]; exit }
    ' "$SAMPLES"
}

patients() { sample_rows | cut -f1 | awk '!seen[$0]++'; }
normals()  { sample_rows | awk -F'\t' '$3 == 0 { print $1 "\t" $2 }'; }   # patient normal
tumours()  { sample_rows | awk -F'\t' '$3 == 1 { print $1 "\t" $2 }'; }   # patient tumour
sex_of()   { sample_rows | awk -F'\t' -v p="$1" '$1 == p { print $4; exit }'; }

# pairs: patient tumour normal -- every tumour against its patient's normal.
# A patient with several normals uses the first one listed (warned).
pairs() {
    sample_rows | awk -F'\t' '
        $3 == 0 { if ($1 in n) dup[$1] = 1; else n[$1] = $2; next }
        $3 == 1 { t[++k] = $1 "\t" $2 }
        END {
            for (p in dup) print "WARNING: " p " has more than one normal; using " n[p] > "/dev/stderr"
            for (i = 1; i <= k; i++) {
                split(t[i], f, "\t")
                if (f[1] in n) print f[1] "\t" f[2] "\t" n[f[1]]
                else print "WARNING: no normal for " f[2] ", skipping" > "/dev/stderr"
            }
        }'
}

# -----------------------------------------------------------------------------
# Where things live. Change a layout here and every stage follows.
# -----------------------------------------------------------------------------
pair_id() { echo "${1}_vs_${2}"; }   # tumour normal -> sarek's naming

recal_cram() {
    local o; o=$(sample_field "$1" cram)
    if [[ -n $o ]]; then echo "$o"; else echo "$SAREK_DIR/preprocessing/recalibrated/$1/$1.recal.cram"; fi
}

# recal_bam: override > bam sarek wrote itself > stage-01 cram_to_bam output
recal_bam() {
    local o; o=$(sample_field "$1" bam)
    if [[ -n $o ]]; then echo "$o"; return; fi
    local s="$SAREK_DIR/preprocessing/recalibrated/$1/$1.recal.bam"
    if [[ -e $s ]]; then echo "$s"; else echo "$BAM_DIR/$1.recal.bam"; fi
}

# mapped_bam: pre-BQSR alignment (sarek --save_mapped); falls back to recal_bam
mapped_bam() {
    local s="$SAREK_DIR/preprocessing/mapped/$1/$1.bam"
    if [[ -e $s ]]; then echo "$s"; else recal_bam "$1"; fi
}

mutect2_raw()      { local id; id=$(pair_id "$1" "$2"); echo "$SAREK_DIR/variant_calling/mutect2/$id/$id.mutect2.vcf.gz"; }
mutect2_filtered() { echo "$MUTECT2_DIR/$(pair_id "$1" "$2").mutect2_filtered.vcf.gz"; }
ascat_cnv()        { local id; id=$(pair_id "$1" "$2"); echo "$SAREK_DIR/variant_calling/ascat/$id/$id.cnvs.txt"; }
ascat_pp()         { local id; id=$(pair_id "$1" "$2"); echo "$SAREK_DIR/variant_calling/ascat/$id/$id.purityploidy.txt"; }

# consensus_vcf <tumour>: stage 03's 2-of-3 consensus VCF for this tumour.
# VCF tumour names can carry suffixes (COS_001_063_PRG_merged), so each VCF
# goes to the LONGEST tumour sample ID its name starts with. Prints nothing
# if there is no match.
consensus_vcf() {
    local v t best
    for v in "$CONSENSUS_DIR"/vcf/*.vcf; do
        [[ -e $v ]] || return 0
        t=$(grep -m1 '^#CHROM' "$v" | cut -f10)
        best=$(tumours | cut -f2 | awk -v t="$t" 'index(t, $0) == 1 && length($0) > length(b) { b = $0 } END { print b }')
        if [[ $best == "$1" ]]; then echo "$v"; return 0; fi
    done
}

# -----------------------------------------------------------------------------
# Job helpers
# -----------------------------------------------------------------------------

# submit <job_name> <n_tasks> <job_script> [args...]
# Always an array job (-t 1-N), so job scripts can rely on SGE_TASK_ID.
# DRY_RUN=1 prints the qsub line instead; HOLD=<jobname> adds -hold_jid.
submit() {
    local name=$1 n=$2 script=$3
    shift 3
    ((n > 0)) || die "nothing to submit for $name (0 tasks)"
    mkdir -p "$LOG_DIR"
    local cmd=(qsub -N "$name" -wd "$LOG_DIR" -t "1-$n"
        -v "PIPELINE_DIR=$PIPELINE_DIR,PIPELINE_CONFIG=$PIPELINE_CONFIG,SAMPLES=$SAMPLES")
    [[ -n ${HOLD:-} ]] && cmd+=(-hold_jid "$HOLD")
    cmd+=("$script" "$@")
    if [[ ${DRY_RUN:-0} == 1 ]]; then
        echo "[dry-run] ${cmd[*]}"
    else
        "${cmd[@]}"
    fi
}

# task_line: the SGE_TASK_ID'th line of stdin (TASK=n overrides, for testing
# a job script interactively: TASK=2 bash stages/05_kir/t1k_kir.job.sh)
task_line() {
    local i=${TASK:-${SGE_TASK_ID:-}}
    [[ $i =~ ^[0-9]+$ ]] || die "no task index (run through run.sh, or set TASK=n)"
    local line; line=$(sed -n "${i}p")
    [[ -n $line ]] || die "task $i is past the end of the list"
    echo "$line"
}

use_conda() {
    # shellcheck disable=SC1090
    source "$CONDA_SH"
    conda activate "$1"
}

# Nextflow head jobs: java + apptainer, with the apptainer pull cache and tmp
# on Scratch. TMPDIR must be set HERE (in the job launching nextflow), not in
# a nextflow config: image pulls happen in the head process before any task
# is submitted (see the "Failure 4" notes in config/nextflow/myriad_pvacseq.config).
nextflow_env() {
    module load "$JAVA_MODULE"
    module load "$APPTAINER_MODULE"
    module load samtools
    export NXF_ENABLE_FS_SYNC=true
    export NXF_SINGULARITY_CACHEDIR=$APPTAINER_CACHE/pull
    export TMPDIR=$APPTAINER_CACHE/tmp
    export APPTAINER_TMPDIR=$TMPDIR SINGULARITY_TMPDIR=$TMPDIR
    mkdir -p "$TMPDIR" "$NXF_SINGULARITY_CACHEDIR"
}

# run_steps <usage> <requested> <step>... : call step_<requested> if it is one
# of the listed steps. Steps are submitted one at a time on purpose: each one
# checks that the previous step's files exist before it queues anything.
run_steps() {
    local usage=$1 want=$2
    shift 2
    if [[ -z $want || $want == -h || $want == --help ]]; then
        echo "$usage" >&2
        exit 1
    fi
    local s
    for s in "$@"; do
        if [[ $s == "$want" ]]; then "step_$s"; return; fi
    done
    die "unknown step '$want'. $usage"
}
