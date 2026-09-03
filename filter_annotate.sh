#!/bin/bash -l
#$ -wd /home/regma03/Scratch/
#$ -l h_rt=5:00:00
#$ -l mem=40G

module load samtools
module load gatk


cd  /myriadfs/home/regma03/Scratch/WES_tumor_batch1/


ref=/home/regma03/Scratch/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
vcf_gnomad=/home/regma03/Scratch/af-only-gnomad.hg38.vcf.gz
pon=/home/regma03/Scratch/somatic-hg38_1000g_pon.hg38.vcf
funcotator_ds=/home/regma03/Scratch/funcotator_dataSources.v1.7.20200521s


gatk FilterMutectCalls \
    -R $ref \
    --min-median-mapping-quality 20 \
    --min-median-base-quality 20 \
    -V /myriadfs/home/regma03/Scratch/WES_tumor_batch1/variant_calling/mutect2/COS_001_063_PRG_vs_COS_001_063_PBMC/COS_001_063_PRG_vs_COS_001_063_PBMC.mutect2.vcf.gz \
    -O COS_001_063_PRG_vs_COS_001_063_PBMC.mutect2_filtered.vcf \
    --min-allele-fraction 0.01

gatk Funcotator \
    -R $ref \
    -V COS_001_063_PRG_vs_COS_001_063_PBMC.mutect2_filtered.vcf \
    -O COS_001_063_PRG_vs_COS_001_063_PBMC.mutect2_filtered.maf \
    --remove-filtered-variants TRUE \
    --three-prime-flank-size 1000 \
    --output-file-format MAF \
    --ref-version hg38 \
    --annotation-default tumor_barcode:COS_001_063_PRG\
    --data-sources-path $funcotator_ds




