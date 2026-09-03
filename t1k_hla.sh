#!/bin/bash -l
#$ -wd /home/regma03/Scratch/
#$ -l h_rt=40:00:00
#$ -l mem=16G
#$ -pe smp 4


cd /myriadfs/home/regma03/Scratch/WES_tumor_batch1/HLA/
source ~/Scratch/miniforge3/etc/profile.d/conda.sh
conda activate t1k
#run-t1k -1 /myriadfs/home/regma03/Scratch/WES_tumor_batch1/COS_001_003_PBMC/COS_001_003_PBMC_EKDN260018615-1A_23HH2TLT4_L3_R1.fastq.gz -2 /myriadfs/home/regma03/Scratch/WES_tumor_batch1/COS_001_003_PBMC/COS_001_003_PBMC_EKDN260018615-1A_23HH2TLT4_L3_R2.fastq.gz -b /myriadfs/home/regma03/Scratch/WES_tumor_batch1/preprocessing/recalibrated/COS_001_003_PBMC/COS_001_003_PBMC.recal.bam --preset kir-wes -f /myriadfs/home/regma03/Scratch/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa -o COS_001_003

#run-t1k -1 /myriadfs/home/regma03/Scratch/WES_tumor_batch1/COS_001_003_PBMC/COS_001_003_PBMC_EKDN260018615-1A_23HH2TLT4_L3_R1.fastq.gz -2 /myriadfs/home/regma03/Scratch/WES_tumor_batch1/COS_001_003_PBMC/COS_001_003_PBMC_EKDN260018615-1A_23HH2TLT4_L3_R2.fastq.gz --preset hla-wgs -f /myriadfs/home/regma03/Scratch/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa -o COS_001_003 -t 4

#run-t1k -b  /myriadfs/home/regma03/Scratch/WES_tumor_batch1/preprocessing/recalibrated/COS_001_003_PBMC/COS_001_003_PBMC.recal.bam --preset hla-wgs -f /myriadfs/home/regma03/Scratch/WES_tumor_batch1/HLA/T1K-1.0.10/hlaidx/hlaidx_dna_seq.fa -c /myriadfs/home/regma03/Scratch/WES_tumor_batch1/HLA/T1K-1.0.10/hlaidx/hlaidx_dna_coord.fa -o COS_001_003_HLA -t 4

run-t1k -b  /myriadfs/home/regma03/Scratch/WES_tumor_batch1/preprocessing/recalibrated/COS_003_009_PBMC/COS_003_009_PBMC.recal.bam --preset hla-wgs -f /myriadfs/home/regma03/Scratch/WES_tumor_batch1/HLA/T1K-1.0.10/kiridx/kiridx_dna_seq.fa -c /myriadfs/home/regma03/Scratch/WES_tumor_batch1/HLA/T1K-1.0.10/kiridx/kiridx_dna_coord.fa -o /myriadfs/home/regma03/Scratch/WES_tumor_batch1/KIR/COS_003_009_KIR -t 4

run-t1k -b  /myriadfs/home/regma03/Scratch/WES_tumor_batch1/preprocessing/recalibrated/COS_001_063_PBMC/COS_001_063_PBMC.recal.bam --preset hla-wgs -f /myriadfs/home/regma03/Scratch/WES_tumor_batch1/HLA/T1K-1.0.10/kiridx/kiridx_dna_seq.fa -c /myriadfs/home/regma03/Scratch/WES_tumor_batch1/HLA/T1K-1.0.10/kiridx/kiridx_dna_coord.fa -o /myriadfs/home/regma03/Scratch/WES_tumor_batch1/KIR/COS_001_063_KIR -t 4

run-t1k -b  /myriadfs/home/regma03/Scratch/WES_tumor_batch1/preprocessing/recalibrated/COS_001_037_PBMC/COS_001_037_PBMC.recal.bam --preset hla-wgs -f /myriadfs/home/regma03/Scratch/WES_tumor_batch1/HLA/T1K-1.0.10/kiridx/kiridx_dna_seq.fa -c /myriadfs/home/regma03/Scratch/WES_tumor_batch1/HLA/T1K-1.0.10/kiridx/kiridx_dna_coord.fa -o /myriadfs/home/regma03/Scratch/WES_tumor_batch1/KIR/COS_001_037_KIR -t 4

run-t1k -b  /myriadfs/home/regma03/Scratch/WES_tumor_batch1/HLA/COS_001_018_G/COS_001_018_G.recal.bam --preset hla-wgs -f /myriadfs/home/regma03/Scratch/WES_tumor_batch1/HLA/T1K-1.0.10/kiridx/kiridx_dna_seq.fa -c /myriadfs/home/regma03/Scratch/WES_tumor_batch1/HLA/T1K-1.0.10/kiridx/kiridx_dna_coord.fa -o /myriadfs/home/regma03/Scratch/WES_tumor_batch1/KIR/COS_001_018_G_KIR -t 4

run-t1k -b  /myriadfs/home/regma03/Scratch/WES_tumor_batch1/HLA/COS_004_003_G/COS_004_003_G.recal.bam --preset hla-wgs -f /myriadfs/home/regma03/Scratch/WES_tumor_batch1/HLA/T1K-1.0.10/kiridx/kiridx_dna_seq.fa -c /myriadfs/home/regma03/Scratch/WES_tumor_batch1/HLA/T1K-1.0.10/kiridx/kiridx_dna_coord.fa -o /myriadfs/home/regma03/Scratch/WES_tumor_batch1/KIR/COS_004_003_KIR -t 4

