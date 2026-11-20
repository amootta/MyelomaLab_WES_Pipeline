# Stage 11 helpers. Sourced (after lib/common.sh) by run.sh and lohhla.job.sh.
# The LOHHLA_* settings live in config/pipeline.config; the fallbacks below
# only keep an older batch config (without them) working.

LOHHLA_DIR=${LOHHLA_DIR:-$BATCH_DIR/LOHHLA}
LOHHLA_HOME=${LOHHLA_HOME:-$REF_DIR/lohhla}
LOHHLA_SCRIPT=${LOHHLA_SCRIPT:-$LOHHLA_HOME/LOHHLAscript.hg38.R}
LOHHLA_HLA_FASTA=${LOHHLA_HLA_FASTA:-$LOHHLA_HOME/imgt/hla_abc_gen.fasta}
LOHHLA_HLA_DAT=${LOHHLA_HLA_DAT:-$LOHHLA_HOME/imgt/hla.dat}
LOHHLA_HLA_SOURCE=${LOHHLA_HLA_SOURCE:-${HLA_SOURCE_FOR_PVACSEQ:-optitype}}
LOHHLA_MIN_COV=${LOHHLA_MIN_COV:-10}
LOHHLA_PLOTS=${LOHHLA_PLOTS:-TRUE}
LOHHLA_CLEANUP=${LOHHLA_CLEANUP:-FALSE}
LOHHLA_CONDA_ENV=${LOHHLA_CONDA_ENV:-lohhla}
LOHHLA_NOVO_MODULE=${LOHHLA_NOVO_MODULE-novocraft/3.04.06}
LOHHLA_NOVO_DIR=${LOHHLA_NOVO_DIR:-$REF_DIR/novocraft}

# mod_load <module>: `module load` that also works from `./wes.sh`. That is
# not a login shell, so `module` is not defined there (jobs use bash -l).
mod_load() {
    if ! type module > /dev/null 2>&1; then
        local f
        for f in "${MODULESHOME:-/nonexistent}/init/bash" /etc/profile.d/modules.sh /usr/share/Modules/init/bash; do
            if [[ -f $f ]]; then source "$f"; break; fi
        done
    fi
    type module > /dev/null 2>&1 || die "no 'module' command -- run with: bash -l $0 ..."
    module load "$1"
}

# novo_setup: load the module (or use LOHHLA_NOVO_DIR), put novoindex on PATH
# (LOHHLA calls it bare) and set novo_dir for --novoDir.
novo_setup() {
    if [[ -n $LOHHLA_NOVO_MODULE ]]; then
        mod_load "$LOHHLA_NOVO_MODULE"
        novo_dir=$(dirname "$(command -v novoalign)")
    else
        novo_dir=$LOHHLA_NOVO_DIR
        export PATH=$novo_dir:$PATH
    fi
    [[ -x $novo_dir/novoalign ]] && command -v novoindex > /dev/null ||
        die "novoalign/novoindex not found (LOHHLA_NOVO_MODULE=$LOHHLA_NOVO_MODULE, LOHHLA_NOVO_DIR=$LOHHLA_NOVO_DIR)"
}
