# MyelomaLab_WES_Pipeline
Collection of tools and pipelines used to analyze paired tumour-normal WES data in the Myeloma Lab

**SAREK PIPELINE**
Files to run: 
run_sarek.sh - run sarek pipeline from fastq (from scratch) or from available bam/cram
samplesheet_fastq.csv - sample sheet needed to run sarek from fastq
samplesheet_cram_sarek.csv - sample sheet needed to run sarek from bam/cram
filter_annotate.sh - script to annotate outputs from mutect2 (from vcfs) using funcotator (GATK pipeline)

Additional files needed:
myriad_extra.config - config file tailored to UCL Myriad (need this!)

Reference files needed:

