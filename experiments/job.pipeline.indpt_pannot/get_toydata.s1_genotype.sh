#!/bin/bash
#SBATCH -c 8
#SBATCH -N 1
#SBATCH -t 0-01:00
#SBATCH --array=1-7
#SBATCH -p short
#SBATCH --mem=64000
#SBATCH -o /home/jz286/WES_analysis/wes_rare/experiments/job_info/hostname_%j.out  
#SBATCH -e /home/jz286/WES_analysis/wes_rare/experiments/job_info/hostname_%j.err 
#SBATCH --mail-type=NONE#SBATCH --mail-type=NONE

INPUT_PATH=/n/groups/price/UKBiobank/WES_50K
OUT_PATH=/n/groups/price/martin/WES_analysis/toy_10K
ID_LIST=/n/groups/price/martin/WES_analysis/toy_10K/ID10K_unrelated.txt
SNP_LIST=/n/groups/price/martin/WES_analysis/toy_10K/SNP50K_CHR1t10_MAFg5_noMHC.txt


for CHROM in {1..10}
# for CHROM in 9
do

# Get PLINK2
plink2 \
    --bfile $INPUT_PATH/chr${CHROM}_v1.SPB.hg19 \
    --keep $ID_LIST \
    --extract $SNP_LIST\
    --geno 0.1 \
    --mac 1 \
    --hwe 1e-50\
    --maj-ref \
    --make-pgen \
    --out ${OUT_PATH}/chr${CHROM}_v1.SPB.hg19.toy_10K\
    && rm ${OUT_PATH}/chr${CHROM}_v1.SPB.hg19.toy_10K.log
    
# Recompute the frq info
plink2 \
    --pfile ${OUT_PATH}/chr${CHROM}_v1.SPB.hg19.toy_10K\
    --freq \
    --out ${OUT_PATH}/chr${CHROM}_v1.SPB.hg19.toy_10K\
    && rm ${OUT_PATH}/chr${CHROM}_v1.SPB.hg19.toy_10K.log
done