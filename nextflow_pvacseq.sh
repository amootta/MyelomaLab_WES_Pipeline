#!/bin/bash -l
#$ -wd /home/regma03/Scratch/
#$ -l h_rt=40:00:00
#$ -l mem=16G
#$ -pe smp 1

module load java/temurin-17/17.0.2_8


module load apptainer
module load samtools
export NXF_ENABLE_FS_SYNC=true
export TMPDIR="$HOME/Scratch/.apptainer/tmp"
export APPTAINER_TMPDIR="$HOME/Scratch/.apptainer/tmp"
export SINGULARITY_TMPDIR="$HOME/Scratch/.apptainer/tmp"
#mkdir -p "$TMPDIR"


#mkdir -p "$HOME/Scratch/.apptainer/pull"

#apptainer pull "$HOME/Scratch/.apptainer/pull/docker.io-jupyter-scipy-notebook-python-3.11.img" docker://docker.io/jupyter/scipy-notebook:python-3.11
#apptainer pull "$HOME/Scratch/.apptainer/pull/docker.io-griffithlab-pvactools-6.0.5.img" docker://docker.io/griffithlab/pvactools:6.0.5
cd /myriadfs/home/regma03/Scratch/mgp1000/output/somatic/consensus/
nextflow run nextflow-pvacseq-v.1.0.0/main.nf --max_time 40h -profile conda --input ./COS_001_063_PRG_vcf/ --hla_csv HLA_pvacseq_COS_001_063.csv --fasta /myriadfs/home/regma03/Scratch/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa --outdir . -c myriad_pvaseq.config --pvacseq_iedb ./iedb/ --vep_plugins ./download/VEP_plugins/ --pvacseq_algorithm 'MHCflurry MHCnuggetsI NetMHC'
