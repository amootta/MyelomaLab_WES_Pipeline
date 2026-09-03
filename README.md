# MyelomaLab_WES_Pipeline
Collection of tools and pipelines used to analyze paired tumour-normal WES data in the Myeloma Lab

Doesn't yet include parallelization..

Other files can be found on the RDSS: 

Please let me know if you need access and I can add you! 


**SAREK PIPELINE**

**Scripts**

run_sarek.sh - run sarek pipeline from fastq (from scratch) or from available bam/cram

filter_annotate.sh - script to annotate outputs from mutect2 (from vcfs) using funcotator (GATK pipeline). Also produces a filtered VCF file (Mutect2 PASS Filter) 

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


**KIR PIPELINE**

KIR Calling is run on T1K or KIRCLE

**Installation**
*T1K*
To install T1K you will need to first install conda : https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

then you will need to create an appropriate conda environment to run T1K: 

conda env create -f t1k_environment.yml 

T1K runs on BAMs, which you already have from previous pipeline runs

*KIRCLE*
To install KIRCLE, download the KIRCLE zip folder that is in RDSS. 

**Scripts**
t1k.sh - script to run T1K for KIR Calling 
KIRCLE.sh  - script to run KIRCLE for KIR Calling 

**Additional files needed**
Download the T1K zip folder on RDSS and thats where the appropriate reference files are. Change the directory in your code to point to there. 

Reference files for KIRCLE are also in the zipped KIRCLE Folder. Change the BLASTDB export path in the script accordingly


**HLA PIPELINE**

HLA Calling is run using T1K, Nextflow HLATyping (Optitype and HLAHD) and HLA-LA

**Installations**

*T1K* - nothing new to install, use the same installation as the one we used to run KIR

*Nextflow HLATyping* - need to download the hlahd.1.7.1.tar.gz from RDSS. Need to download the config files to run. 

*HLA-LA*  - need to use conda again to install this:

conda env create -f hlala_environment.yml 

AFTERWARDS (slightly complicated). After having installed the conda environment locally, please copy this folder into the appropriate path as below (change to your user/where the environment is installed locally on your myriad user login):

/home/regma03/Scratch/conda/envs/hla-la/opt/hla-la/graphs/PRG_MHC_GRCh38_withIMGT

The PRG_MHC_GRCh38_withIMGT folder is in RDSS for you to download and add to that exact path


**Scripts**
t1k_hla.sh - script to run T1K for HLA calling

nextflow_HLA.sh - script to run the Nextflow HLAtyping pipeline using Optitype and HLAHD

hla_la.sh - script to run HLA-LA pipeline 


**Additional files needed:**

hlatyping_fix.config - config file tailored to UCL Myriad (need this!)for Nextflow HLATyping pipeline

samplesheet_hla.csv - sample sheet needed to run hla calling from BAMs for Nextflow HLATyping pipeline


##########################################################################################################################

That is basically the end of the basic calling pipelines... now onto some downstream stuff! 

**Pyclone**
Pyclone is a tool used to infer the clonal population structure of heterogeneous cancer samples using deep next-generation sequencing (NGS) data

Questions: to use **filtered** somatic mutation calls based on the MGP1000 pipeline where we only use mutations that are called in 2/3 tools - between Strelka, Mutect2 and Varscan? 

Or just use all filtered mutations from Mutect2 only - based on GATK Gold Standards??????

**Installation**

Again need to use conda to install an environment:

**Scripts**
Input needed is a **vcf file of somatic mutations, purityploidy file from ASCAT output, and also cnv file from ASCAT output (these were run using Sarek)**. Need to run this on both BSL and Followup/PRG timepoint together. Make sure you have all of these available.

1. First run pyclone_counts.sh
2. Run pileup.sh on the resulting files from step 1
3. Run  build_pyclone_input.sh with Rscript build_pyclone_input.R
4. Run pyclone_run.sh on the output from 4
5. Summarize and visualize plots of the results from 4 using summarize_pyclone.sh with Rscript summarize_pyclone.R 





**Mutational Signature Calling**

Again the same question of what set of mutations to use: to use **filtered** somatic mutation calls based on the MGP1000 pipeline where we only use mutations that are called in 2/3 tools - between Strelka, Mutect2 and Varscan? 

Or just use all filtered mutations from Mutect2 only - based on GATK Gold Standards??????








