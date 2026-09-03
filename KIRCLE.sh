#!/bin/bash -l
#$ -wd /myriadfs/home/regma03/Scratch/KIRCLE/KIRCLE
#$ -l h_rt=40:00:00
#$ -l mem=16G
#$ -pe smp 1

cd /myriadfs/home/regma03/Scratch/KIRCLE/KIRCLE
module unload gcc-libs
module load blast+
module load python3
export BLASTDB=/myriadfs/home/regma03/Scratch/KIRCLE/ref_files/ref_fastas/
module load samtools

python /myriadfs/home/regma03/Scratch/KIRCLE/KIRCLE/KIRCLE.py -i /myriadfs/home/regma03/Scratch/WES_tumor_batch1/preprocessing/recalibrated/COS_001_003_PBMC/COS_001_003_PBMC.recal.bam -o /myriadfs/home/regma03/Scratch/WES_tumor_batch1/KIR/COS_001_003 -l /myriadfs/home/regma03/Scratch/KIRCLE/ref_files/hg38_KIR_locations_noKV_noalt.tsv -g hg38