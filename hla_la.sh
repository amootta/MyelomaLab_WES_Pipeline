#!/bin/bash -l
#$ -wd /home/regma03/Scratch/
#$ -l h_rt=40:00:00
#$ -l mem=16G
#$ -pe smp 4

cd /myriadfs/home/regma03/Scratch/WES_tumor_batch1/HLA/
source ~/Scratch/miniforge3/etc/profile.d/conda.sh

conda activate hla-la
#HLA-LA.pl --BAM /myriadfs/home/regma03/Scratch/WES_tumor_batch1/HLA/COS_004_003_G/COS_004_003_G.recal.bam --graph PRG_MHC_GRCh38_withIMGT --sampleID COS_004_003 --maxThreads 4 --workingDir /myriadfs/home/regma03/Scratch/WES_tumor_batch1/HLA/

HLA-LA.pl --BAM /myriadfs/home/regma03/Scratch/WES_tumor_batch1/preprocessing/recalibrated/COS_003_009_PBMC/COS_003_009_PBMC.recal.bam --graph PRG_MHC_GRCh38_withIMGT --sampleID COS_003_009 --maxThreads 4 --workingDir /myriadfs/home/regma03/Scratch/WES_tumor_batch1/HLA/HLA_LA

HLA-LA.pl --BAM /myriadfs/home/regma03/Scratch/WES_tumor_batch1/preprocessing/recalibrated/COS_001_063_PBMC/COS_001_063_PBMC.recal.bam --graph PRG_MHC_GRCh38_withIMGT --sampleID COS_001_063 --maxThreads 4 --workingDir /myriadfs/home/regma03/Scratch/WES_tumor_batch1/HLA/HLA_LA

HLA-LA.pl --BAM /myriadfs/home/regma03/Scratch/WES_tumor_batch1/preprocessing/recalibrated/COS_001_037_PBMC/COS_001_037_PBMC.recal.bam --graph PRG_MHC_GRCh38_withIMGT --sampleID COS_001_037 --maxThreads 4 --workingDir /myriadfs/home/regma03/Scratch/WES_tumor_batch1/HLA/HLA_LA

HLA-LA.pl --BAM /myriadfs/home/regma03/Scratch/WES_tumor_batch1/HLA/COS_001_018_G/COS_001_018_G.recal.bam --graph PRG_MHC_GRCh38_withIMGT --sampleID COS_001_018 --maxThreads 4 --workingDir /myriadfs/home/regma03/Scratch/WES_tumor_batch1/HLA/HLA_LA

HLA-LA.pl --BAM /myriadfs/home/regma03/Scratch/WES_tumor_batch1/preprocessing/recalibrated/COS_001_003_PBMC/COS_001_003_PBMC.recal.bam --graph PRG_MHC_GRCh38_withIMGT --sampleID COS_001_003 --maxThreads 4 --workingDir /myriadfs/home/regma03/Scratch/WES_tumor_batch1/HLA/HLA_LA
