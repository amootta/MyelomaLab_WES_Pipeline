# MyelomaLab_WES_Pipeline
Collection of tools and pipelines used to analyze paired tumour-normal WES data in the Myeloma Lab

Doesn't yet include parallelization..

**SAREK PIPELINE**

**Scripts**

run_sarek.sh - run sarek pipeline from fastq (from scratch) or from available bam/cram

filter_annotate.sh - script to annotate outputs from mutect2 (from vcfs) using funcotator (GATK pipeline)

**Additional files needed:**

myriad_extra.config - config file tailored to UCL Myriad (need this!)

samplesheet_fastq.csv - sample sheet needed to run sarek from fastq

samplesheet_cram_sarek.csv - sample sheet needed to run sarek from bam/cram

**Reference files needed:**
HG38 reference genome to use is uploaded on to RDSS
Other reference files needed are provided in References folder in RDSS
for Funcotator: download funcotator data source folder: https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/funcotator?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))

Also available on RDSS

**MGP1000 PIPELINE**

Two choices: can either run from the .bams that have been produced by Sarek pipeline, or again start from fastqs (from Scratch)

Personally I prefer starting from the crams/bams that have been produced by Sarek as they follow the GATK Gold Standard pre-processing pipeline, but again personal choice! Have given both versions of the scripts in the files uploaded

MGP1000 installation folder to use that is appropriate with UCL Myriad system has been uploaded onto here as well. Again needs an installation of Nextflow before you run this. Particularly I used Nextflow version 22.10.8 

**Scripts**
mgp1000_WES.sh - full run script for MGP1000. If using .crams produced by Sarek, includes converting .crams to .bams as MGP1000 only takes .bams


**Reference files needed:**

Other reference files needed are provided in References folder in RDSS





