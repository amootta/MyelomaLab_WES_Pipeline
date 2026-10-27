#!/bin/bash -l
#$ -wd /home/regma03/Scratch/
#$ -l h_rt=5:00:00

module load samtools
module load gatk

cd /home/regma03/Scratch/CARDAMON_mutect2


ref=/home/regma03/Scratch/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
vcf_gnomad=/home/regma03/Scratch/af-only-gnomad.hg38.vcf.gz
pon=/home/regma03/Scratch/somatic-hg38_1000g_pon.hg38.vcf
funcotator_ds=/home/regma03/Scratch/funcotator_dataSources.v1.7.20200521s

# get VCF corresponding to array index
VCF=$1

# remove extension to create prefix
OUT_PREFIX=$(basename "$VCF" .vcf)

# optional: extract sample ID
SAMPLEID=${OUT_PREFIX%%_*}


echo "Processing $VCF"
echo "Sample: $SAMPLEID"


gatk FilterMutectCalls \
    -R $ref \
    --min-median-mapping-quality 20 \
    --min-median-base-quality 20 \
    -V "$VCF" \
    -O "${OUT_PREFIX}_filtered.vcf" \
    --min-allele-fraction 0.01

#gatk Funcotator \
#    -R $ref \
#    -V "${OUT_PREFIX}_filtered.vcf" \
#    -O "${OUT_PREFIX}_filtered.maf" \
#    --remove-filtered-variants TRUE \
#    --three-prime-flank-size 1000 \
#    --output-file-format MAF \
#    --ref-version hg38 \
#    --annotation-default tumor_barcode:"${SAMPLEID}_T" \
#    --data-sources-path $funcotator_ds