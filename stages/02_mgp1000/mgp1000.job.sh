#!/bin/bash -l
#$ -l h_rt=48:00:00
#$ -l mem=16G
#$ -pe smp 2
# Nextflow head job for one mgp1000 workflow (it submits its own SGE jobs via
# --executor sge). From mgp1000_WES.sh. The head job only orchestrates, but it
# has to outlive the whole run; if it hits the wall-time limit, rerun the same
# wes.sh command and -resume carries on.
# usage: mgp1000.job.sh <somatic|germline> <sample sheet>
source "$PIPELINE_DIR/lib/common.sh"
workflow=$1 sheet=$2

for m in $MGP1000_APPTAINER_MODULES; do module load "$m"; done
module load samtools
module load "$JAVA_MODULE"
export NXF_VER=$MGP1000_NXF_VER
export NXF_SINGULARITY_CACHEDIR=$MGP1000_SINGULARITY_CACHE

case $workflow in
    somatic)  tools=$MGP1000_SOMATIC_TOOLS ;;
    germline) tools=$MGP1000_GERMLINE_TOOLS ;;
    *) die "unknown mgp1000 workflow '$workflow'" ;;
esac

# The install in MGP1000_DIR is shared by every batch; launching from this
# batch's MGP1000_RUN_DIR keeps work/, .nextflow/ and output/ per batch.
mkdir -p "$MGP1000_RUN_DIR"
cd "$MGP1000_RUN_DIR"
# mgp1000 names its images and reference files by relative path
# (containers/<tool>.simg, references/hg38/...), which Nextflow looks for in
# the launch dir, so link the install's folders here.
[[ -e containers ]] || ln -s "$MGP1000_DIR/containers" containers
# references/ is this batch's own copy of links into the install, except that
# Homo_sapiens_assembly38.fasta/.fai/.dict point at REF_FASTA: the bams come
# from sarek on REF_FASTA, and Manta/Strelka/GATK refuse a FASTA with contigs
# the bams lack (the GATK hg38 has ~3,400, the no-alt set 195). The bwa index
# (.fasta.64.*) is left out; only SvABA uses it, and it is off.
[[ -L references ]] && rm references   # older runs linked the whole folder
for d in "$MGP1000_DIR"/references/*/; do
    d=${d%/} sub=references/$(basename "$d")
    mkdir -p "$sub"
    for f in "$d"/*; do
        case $(basename "$f") in Homo_sapiens_assembly38.fasta*|Homo_sapiens_assembly38.dict) continue ;; esac
        ln -sfn "$f" "$sub/$(basename "$f")"
    done
done
for f in "$MGP1000_DIR"/references/*; do
    [[ -d $f ]] || ln -sfn "$f" "references/$(basename "$f")"
done
need "$REF_FASTA" "$REF_FAI"
ln -sfn "$REF_FASTA" references/hg38/Homo_sapiens_assembly38.fasta
ln -sfn "$REF_FAI" references/hg38/Homo_sapiens_assembly38.fasta.fai
dict=references/hg38/Homo_sapiens_assembly38.dict
if [[ -e ${REF_FASTA%.*}.dict ]]; then
    ln -sfn "${REF_FASTA%.*}.dict" "$dict"
elif [[ ! -f $dict || -L $dict ]]; then   # no dict beside REF_FASTA: make one here
    rm -f "$dict"
    samtools dict -o "$dict" "$REF_FASTA"
fi
# sarek names each read group's sample <patient>_<sample>, but mgp1000 passes
# the bare sample (e.g. Mutect2 --normal-sample), so give each staged bam that
# needs it a new header and index. These are new files that replace the hard
# links, so sarek's own bam and .bai are never written through. somatic and
# germline may start together, hence the lock.
exec 9> "$MGP1000_INPUT_DIR/.rename_sm.lock"
flock 9
for bam in "$MGP1000_INPUT_DIR"/*.final.bam; do
    s=$(basename "$bam" .final.bam)
    samtools view -H "$bam" | awk -F'\t' -v sm="SM:$s" '
        /^@RG/ { for (i = 2; i <= NF; i++) if ($i ~ /^SM:/ && $i != sm) bad = 1 }
        END { exit !bad }' || continue
    log "$s: setting read-group SM to $s in $bam"
    samtools view -H "$bam" | awk -v sm="SM:$s" 'BEGIN { FS = OFS = "\t" }
        /^@RG/ { for (i = 2; i <= NF; i++) if ($i ~ /^SM:/) $i = sm } 1' > "$bam.hdr"
    samtools reheader "$bam.hdr" "$bam" > "$bam.tmp"
    mv "$bam.tmp" "$bam"
    rm -f "$bam.hdr" "$bam.bai"
    samtools index -@ 2 "$bam"
done
flock -u 9
# somatic and germline share this launch dir, so a bare -resume would pick
# whichever ran last (and fail on its lock if it is still running). Resume
# this workflow's own last session instead (history: ... session id, command).
# A session belongs to the workflow it was first started with: a germline run
# that resumed somatic's session by mistake is logged under somatic's id.
sid=$(awk -F'\t' -v nf="/$workflow.nf" '!($6 in first) { first[$6] = $7 } index(first[$6], nf) { s = $6 } END { print s }' .nextflow/history 2>/dev/null)
resume=(-resume ${sid:+"$sid"})
[[ -z $sid && -s .nextflow/history ]] && resume=()   # first run of this workflow here
# shellcheck disable=SC2086  # $tools is a list of flags
nextflow run "$MGP1000_DIR/$workflow.nf" \
    --run_id "${MGP1000_RUN_ID}_$workflow" \
    --seq_protocol WES \
    --sample_sheet "$sheet" \
    --input_dir "$MGP1000_INPUT_DIR" \
    --output_dir "$MGP1000_RUN_DIR/output" \
    --email "$EMAIL" \
    --executor sge \
    -profile "$workflow" \
    $tools \
    "${resume[@]}"

# --- not wired into wes.sh ---------------------------------------------------
# Preprocessing from FASTQs (instead of sarek's bams), once per sample, before
# germline/somatic. FASTQs must be named *_R1/_R2 and end in .fq.gz:
#   nextflow run preprocessing.nf --run_id <sample> --input_format fastq \
#     --input_dir <fastq dir> --email "$EMAIL" -profile preprocessing \
#     --executor sge --seq_protocol WES --cpus 8 --output_dir .
