#!/bin/bash -l
#$ -l h_rt=4:00:00
#$ -l mem=8G
#$ -pe smp 1
# Per-caller mgp1000 SNV VCFs -> bcftools norm -> per-caller hard filters ->
# sites called by >= 2 of the 3 callers -> per-pair variant table for MMSig.
# Port of merging_results.sh (thanks Filippo), with two fixes:
#   * file names are normalised to <pair>.<caller>.*, so all three callers line
#     up (the old strelka glob never matched its own normalised file names);
#   * the 2-of-3 set is taken from `bcftools isec -n +2` sites directly. The old
#     `-w 1` wrote only Mutect2's records, silently dropping strelka+varscan-only sites.
source "$PIPELINE_DIR/lib/common.sh"
module load samtools
module load bcftools

out=$CONSENSUS_DIR/isec
norm=$out/normalized
mkdir -p "$norm" "$out/variant_tsv"
log_tn=$out/sample_order_check.log
echo -e "vcf\tcaller\ttumor_sample\ttumor_idx\tnormal_sample\tnormal_idx" > "$log_tn"

# <pair> from e.g. T_vs_N.mutect.somatic.snv.vcf.gz / T_vs_N.somatic.snv.vcf.gz
# TODO: check against real mgp1000 file names -- this assumes the pair ID is
#       everything before the first of .mutect/.strelka/.varscan/.somatic.
pair_of() { basename "$1" .vcf.gz | sed -E 's/\.(mutect|strelka|varscan|somatic|snv)(\..*)?$//'; }

# ---- normalise --------------------------------------------------------------
for c in mutect strelka varscan; do
    for vcf in "$MGP1000_SOMATIC_DIR/$c"/*.somatic.snv.vcf.gz; do
        [[ -e $vcf ]] || { warn "no $c VCFs"; break; }
        p=$(pair_of "$vcf")
        bcftools norm -m -any -f "$REF_FASTA" "$vcf" -Oz -o "$norm/$p.$c.norm.vcf.gz"
        bcftools index -f "$norm/$p.$c.norm.vcf.gz"
    done
done

# get_tn_indices <vcf> <caller> -> "tumor_idx tumor_name normal_idx normal_name"
# Mutect2's ##tumor_sample/##normal_sample headers are ground truth; otherwise
# fall back to TUMOR/NORMAL-style names. Never guesses silently.
get_tn_indices() {
    local vcf=$1 caller=$2 samples tumor_name="" normal_name="" s i
    mapfile -t samples < <(bcftools query -l "$vcf")
    if [[ ${#samples[@]} -ne 2 ]]; then
        echo "ERROR: $vcf has ${#samples[@]} samples, expected 2" >&2
        return 1
    fi
    if [[ $caller == mutect ]]; then
        tumor_name=$(bcftools view -h "$vcf" | grep -m1 '^##tumor_sample=' | sed 's/^##tumor_sample=//')
        normal_name=$(bcftools view -h "$vcf" | grep -m1 '^##normal_sample=' | sed 's/^##normal_sample=//')
    fi
    if [[ -z $tumor_name || -z $normal_name ]]; then
        tumor_name="" normal_name=""
        for s in "${samples[@]}"; do
            case "$s" in
                TUMOR|*_T|*-T|*[Tt]umou?r*) tumor_name=$s ;;
                NORMAL|*_N|*-N|*[Nn]ormal*) normal_name=$s ;;
            esac
        done
    fi
    if [[ -z $tumor_name || -z $normal_name || $tumor_name == "$normal_name" ]]; then
        echo "ERROR: could not resolve tumor/normal for $vcf (samples: ${samples[*]})" >&2
        return 1
    fi
    local t=-1 n=-1
    for i in "${!samples[@]}"; do
        [[ ${samples[$i]} == "$tumor_name" ]] && t=$i
        [[ ${samples[$i]} == "$normal_name" ]] && n=$i
    done
    echo "$t $tumor_name $n $normal_name"
}

# ---- hard filters (thresholds unchanged from merging_results.sh) --------------
declare -A expr
expr[mutect]='FORMAT/AD[@T:1]>=6 && FORMAT/DP[@T]>=30 && INFO/MMQ>=40'
expr[strelka]='((ALT="A" && FORMAT/AU[@T:0]>=6) || (ALT="C" && FORMAT/CU[@T:0]>=6) || (ALT="G" && FORMAT/GU[@T:0]>=6) || (ALT="T" && FORMAT/TU[@T:0]>=6)) && FORMAT/DP[@T]>=30 && INFO/MQ>=40'
expr[varscan]='FORMAT/DP[@T]>=30 && FORMAT/AD[@T]>=5'

for c in mutect strelka varscan; do
    for vcf in "$norm"/*."$c".norm.vcf.gz; do
        [[ -e $vcf ]] || break
        res=$(get_tn_indices "$vcf" "$c") || { warn "skipping $vcf"; continue; }
        read -r t_idx t_name n_idx n_name <<< "$res"
        echo -e "$(basename "$vcf")\t$c\t$t_name\t$t_idx\t$n_name\t$n_idx" >> "$log_tn"
        bcftools filter -i "${expr[$c]//@T/$t_idx}" "$vcf" -Oz -o "${vcf%.norm.vcf.gz}.filtered.vcf.gz"
        bcftools index -f "${vcf%.norm.vcf.gz}.filtered.vcf.gz"
    done
done
log "tumour/normal detection log: $log_tn -- check the varscan rows before trusting the output"

# ---- >= 2 of 3 --------------------------------------------------------------
for m in "$norm"/*.mutect.filtered.vcf.gz; do
    [[ -e $m ]] || die "no filtered mutect VCFs in $norm"
    p=$(basename "$m" .mutect.filtered.vcf.gz)
    s=$norm/$p.strelka.filtered.vcf.gz v=$norm/$p.varscan.filtered.vcf.gz
    [[ -e $s && -e $v ]] || { warn "$p: missing strelka or varscan, skipped"; continue; }
    {
        echo -e "SAMPLE\tCHROM\tPOS\tREF\tALT\tcallers_mutect_strelka_varscan"
        bcftools isec -n +2 "$m" "$s" "$v" | awk -v S="$p" 'BEGIN { OFS = "\t" } { print S, $1, $2, $3, $4, $5 }'
    } > "$out/variant_tsv/$p.snv.filtered.variants.tsv"
    log "$p: $(($(wc -l < "$out/variant_tsv/$p.snv.filtered.variants.tsv") - 1)) consensus SNVs"
done
