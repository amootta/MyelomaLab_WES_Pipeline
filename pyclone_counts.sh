REF=/home/regma03/Scratch/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa


bcftools norm -m -any -f "$REF" "COS_001_063_PBMC_BL_vs_COS_001_063_PBMC.mutect2_filtered.vcf" | bcftools view -f PASS -v snps | bcftools query -s"COS_001_063_COS_001_063_PBMC_BL" --force-samples -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD{0}\t%AD{1}]\n' > ./COS_001_063_BSL.tsv

bcftools norm -m -any -f "$REF" "COS_001_063_PRG_vs_COS_001_063_PBMC.mutect2_filtered.vcf" | bcftools view -f PASS -v snps | bcftools query -s"COS_001_063_COS_001_063_PBMC_PRG" --force-samples -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD{0}\t%AD{1}]\n' > ./COS_001_063_PRG.tsv