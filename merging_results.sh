#!/bin/bash -l
module load samtools
module load bcftools
set -uo pipefail

SOMATIC_DIR=/myriadfs/home/regma03/Scratch/mgp1000/output/somatic
REF=/myriadfs/home/regma03/Scratch/mgp1000/references/hg38/Homo_sapiens_assembly38.fasta
NORM_DIR="${SOMATIC_DIR}/normalized"
mkdir -p "$NORM_DIR"

#############################################
# NORMALIZATION (unchanged from your script)
#############################################
for vcf in "${SOMATIC_DIR}"/strelka/*.somatic.snv.vcf.gz; do
    base=$(basename "$vcf" .vcf.gz)
    bcftools norm -m -any -f "$REF" "$vcf" -Oz \
        -o "${NORM_DIR}/${base}.strelka.norm.vcf.gz"
    bcftools index "${NORM_DIR}/${base}.strelka.norm.vcf.gz"
done

for vcf in "${SOMATIC_DIR}"/mutect/*.mutect.somatic.snv.vcf.gz; do
    base=$(basename "$vcf" .vcf.gz)
    bcftools norm -m -any -f "$REF" "$vcf" -Oz \
        -o "${NORM_DIR}/${base}.mutect.norm.vcf.gz"
    bcftools index "${NORM_DIR}/${base}.mutect.norm.vcf.gz"
done

for vcf in "${SOMATIC_DIR}"/varscan/*.somatic.snv.vcf.gz; do
    base=$(basename "$vcf" .vcf.gz)
    bcftools norm -m -any -f "$REF" "$vcf" -Oz \
        -o "${NORM_DIR}/${base}.varscan.norm.vcf.gz"
    bcftools index "${NORM_DIR}/${base}.varscan.norm.vcf.gz"
done

#############################################
# HARD FILTERING - with tumor/normal auto-detection
#############################################
cd "$NORM_DIR"

LOGFILE="${NORM_DIR}/sample_order_check.log"
echo -e "vcf\tcaller\ttumor_sample\ttumor_idx\tnormal_sample\tnormal_idx" > "$LOGFILE"

# get_tn_indices <vcf> <caller>
# Prints "tumor_idx tumor_name normal_idx normal_name" to stdout on success.
# Returns 1 (prints nothing to stdout) if it can't unambiguously resolve both
# roles - never guesses silently.
#
# Detection strategy, in order of trust:
#   1. mutect: read the ##tumor_sample= / ##normal_sample= header lines that
#      GATK writes itself when Mutect2 is run with -normal <name> (paired
#      mode) - this is ground truth, independent of column order.
#   2. Fallback (also used for strelka/varscan): match sample names against
#      the literal "TUMOR"/"NORMAL" convention Strelka and VarScan default
#      to, or common suffix patterns (_T/_N, -T/-N, "tumor(ur)"/"normal").
get_tn_indices () {
    local vcf="$1"
    local caller="$2"
    local samples
    mapfile -t samples < <(bcftools query -l "$vcf")

    if [[ ${#samples[@]} -ne 2 ]]; then
        echo "ERROR: $vcf does not have exactly 2 samples (${#samples[*]:-0} found: ${samples[*]:-none})" >&2
        return 1
    fi

    local tumor_name="" normal_name=""

    if [[ "$caller" == "mutect" ]]; then
        tumor_name=$(bcftools view -h "$vcf" | grep -m1 '^##tumor_sample=' | sed 's/^##tumor_sample=//')
        normal_name=$(bcftools view -h "$vcf" | grep -m1 '^##normal_sample=' | sed 's/^##normal_sample=//')
    fi

    if [[ -z "$tumor_name" || -z "$normal_name" ]]; then
        tumor_name=""; normal_name=""
        for s in "${samples[@]}"; do
            case "$s" in
                TUMOR|*_T|*-T|*[Tt]umou?r*) tumor_name="$s" ;;
                NORMAL|*_N|*-N|*[Nn]ormal*) normal_name="$s" ;;
            esac
        done
    fi

    if [[ -z "$tumor_name" || -z "$normal_name" || "$tumor_name" == "$normal_name" ]]; then
        echo "ERROR: could not unambiguously resolve tumor/normal for $vcf (samples: ${samples[*]})" >&2
        return 1
    fi

    local tumor_idx=-1 normal_idx=-1 i
    for i in "${!samples[@]}"; do
        [[ "${samples[$i]}" == "$tumor_name" ]] && tumor_idx=$i
        [[ "${samples[$i]}" == "$normal_name" ]] && normal_idx=$i
    done

    if [[ $tumor_idx -eq -1 || $normal_idx -eq -1 ]]; then
        echo "ERROR: resolved names ($tumor_name / $normal_name) don't match sample list for $vcf" >&2
        return 1
    fi

    echo "$tumor_idx $tumor_name $normal_idx $normal_name"
}

# ---- Mutect2 ----
for mutect_vcf in *mutect.somatic*norm.vcf.gz; do
    result=$(get_tn_indices "$mutect_vcf" mutect) || { echo "SKIPPING $mutect_vcf (see error above)"; continue; }
    read -r t_idx t_name n_idx n_name <<< "$result"
    echo -e "${mutect_vcf}\tmutect\t${t_name}\t${t_idx}\t${n_name}\t${n_idx}" >> "$LOGFILE"

    bcftools filter \
        -i "FORMAT/AD[${t_idx}:1]>=6 && FORMAT/DP[${t_idx}]>=30 && INFO/MMQ>=40" \
        "$mutect_vcf" \
        -Oz \
        -o "${mutect_vcf%.norm.vcf.gz}.filtered.vcf.gz"
done

# ---- Strelka2 ----
for strelka_vcf in *strelka.somatic.snv.norm.vcf.gz; do
    result=$(get_tn_indices "$strelka_vcf" strelka) || { echo "SKIPPING $strelka_vcf (see error above)"; continue; }
    read -r t_idx t_name n_idx n_name <<< "$result"
    echo -e "${strelka_vcf}\tstrelka\t${t_name}\t${t_idx}\t${n_name}\t${n_idx}" >> "$LOGFILE"

    bcftools filter \
        -i "((ALT=\"A\" && FORMAT/AU[${t_idx}:0]>=6) || (ALT=\"C\" && FORMAT/CU[${t_idx}:0]>=6) || (ALT=\"G\" && FORMAT/GU[${t_idx}:0]>=6) || (ALT=\"T\" && FORMAT/TU[${t_idx}:0]>=6)) && FORMAT/DP[${t_idx}]>=30 && INFO/MQ>=40" \
        "$strelka_vcf" \
        -Oz \
        -o "${strelka_vcf%.norm.vcf.gz}.filtered.vcf.gz"
done

# ---- VarScan ----
for varscan_vcf in *varscan.somatic.snv.varscan.norm.vcf.gz; do
    result=$(get_tn_indices "$varscan_vcf" varscan) || { echo "SKIPPING $varscan_vcf (see error above)"; continue; }
    read -r t_idx t_name n_idx n_name <<< "$result"
    echo -e "${varscan_vcf}\tvarscan\t${t_name}\t${t_idx}\t${n_name}\t${n_idx}" >> "$LOGFILE"

    bcftools filter \
        -i "FORMAT/DP[${t_idx}]>=30 && FORMAT/AD[${t_idx}]>=5" \
        "$varscan_vcf" \
        -Oz \
        -o "${varscan_vcf%.norm.vcf.gz}.filtered.vcf.gz"
done

echo "Tumor/normal detection log written to: $LOGFILE"
echo "Check it before trusting the filtered output - especially the varscan rows,"
echo "since VarScan sample naming can't be verified against a header tag the way Mutect2's can."

#############################################
# INDEX AFTER FILTERING
#############################################
for vcf in *filtered.vcf.gz; do
    bcftools index "$vcf"
done

#############################################
# CONSENSUS VARIANTS: detected by >=2 of 3 callers (unchanged)
#############################################
mkdir -p consensus variant.tsv
for mutect_vcf in *mutect.somatic.snv.mutect.filtered.vcf.gz; do
    sample=$(basename "$mutect_vcf" .mutect.somatic.snv.mutect.filtered.vcf.gz)
    echo "Processing: $sample"
    strelka_vcf="${sample}.strelka.somatic.snv.filtered.vcf.gz"
    varscan_vcf="${sample}.varscan.somatic.snv.varscan.filtered.vcf.gz"
    if [[ ! -f "$strelka_vcf" ]]; then
        echo "ERROR: Missing Strelka2 file for $sample"
        continue
    fi
    if [[ ! -f "$varscan_vcf" ]]; then
        echo "ERROR: Missing VarScan file for $sample"
        continue
    fi
    bcftools isec \
        -n +2 \
        -w 1 \
        -Oz \
        -o "consensus/${sample}.somatic.snv.consensus.vcf.gz" \
        "$mutect_vcf" \
        "$strelka_vcf" \
        "$varscan_vcf"
    bcftools index -t \
        "consensus/${sample}.somatic.snv.consensus.vcf.gz"
done

for file in consensus/*.somatic.snv.consensus.vcf.gz; do
    sample=$(basename "$file" .somatic.snv.consensus.vcf.gz)
    echo "Creating variant table for: $sample"
    {
        echo -e "SAMPLE\tCHROM\tPOS\tREF\tALT"
        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\n' \
            "$file" |
        awk -v S="$sample" \
            'BEGIN {OFS="\t"} {print S,$0}'
    } > "variant.tsv/${sample}.snv.filtered.variants.tsv"
done