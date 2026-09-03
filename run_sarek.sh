#!/bin/bash -l
#$ -wd /home/regma03/Scratch/
#$ -l h_rt=40:00:00
#$ -l mem=16G
#$ -pe smp 1

########SCRIPT TO RUN SAREK PIPELINE ####################
## Please change directories accordingly
## Uploaded all necessary reference files as well 
## Need to install nextflow manually: https://docs.seqera.io/nextflow/install

module load java/temurin-17/17.0.2_8

module load samtools
export NXF_ENABLE_FS_SYNC=true

cd /myriadfs/home/regma03/Scratch/WES_tumor_batch1/


### if running from fastqs
nextflow run nf-core/sarek -profile ucl_myriad --input samplesheet_cram_sarek.csv --outdir ./ --wes --trim_fastq --igenomes_ignore --fasta /myriadfs/home/regma03/Scratch/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa --fasta_fai /myriadfs/home/regma03/Scratch/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai --dbsnp /myriadfs/home/regma03/Scratch/funcotator_dataSources.v1.8.hg38.20230908s/dbsnp/hg38/hg38_All_20180418.vcf.gz --known_indels /myriadfs/home/regma03/Scratch/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --tools mutect2,ascat,cnvkit,strelka,manta --step variant_calling --email regma03@ucl.ac.uk --save_mapped --ascat_genome hg38 --ascat_loci /myriadfs/home/regma03/Scratch/battenberg_reference/1000G_loci_hg38/battenberg_loci_on_target_hg38/battenberg_loci_on_target_hg38.zip --ascat_alleles /myriadfs/home/regma03/Scratch/battenberg_reference/1000G_loci_hg38/battenberg_alleles_on_target_hg38/battenberg_alleles_on_target_hg38.zip --ascat_loci_gc /home/regma03/Scratch/GC_G1000_on_target_hg38.txt.zip --ascat_loci_rt /home/regma03/Scratch/RT_G1000_on_target_hg38.txt.zip --mappability /home/regma03/Scratch/out100m2_hg38.gem --chr_dir /home/regma03/Scratch/chrfiles_controlfreec -c myriad_extra.config --max_cpus 8

### if running from cram or bams
nextflow run nf-core/sarek -profile ucl_myriad --input samplesheet_fastq.csv --outdir ./ --wes --trim_fastq --igenomes_ignore --fasta /myriadfs/home/regma03/Scratch/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa --fasta_fai /myriadfs/home/regma03/Scratch/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai --dbsnp /myriadfs/home/regma03/Scratch/funcotator_dataSources.v1.8.hg38.20230908s/dbsnp/hg38/hg38_All_20180418.vcf.gz --known_indels /myriadfs/home/regma03/Scratch/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --tools mutect2,ascat,cnvkit,strelka,manta  --email regma03@ucl.ac.uk --save_mapped --ascat_genome hg38 --ascat_loci /myriadfs/home/regma03/Scratch/battenberg_reference/1000G_loci_hg38/battenberg_loci_on_target_hg38/battenberg_loci_on_target_hg38.zip --ascat_alleles /myriadfs/home/regma03/Scratch/battenberg_reference/1000G_loci_hg38/battenberg_alleles_on_target_hg38/battenberg_alleles_on_target_hg38.zip --ascat_loci_gc /home/regma03/Scratch/GC_G1000_on_target_hg38.txt.zip --ascat_loci_rt /home/regma03/Scratch/RT_G1000_on_target_hg38.txt.zip --mappability /home/regma03/Scratch/out100m2_hg38.gem --chr_dir /home/regma03/Scratch/chrfiles_controlfreec -c myriad_extra.config --max_cpus 8
