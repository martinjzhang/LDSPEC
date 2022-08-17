#!/bin/bash
#SBATCH -c 8
#SBATCH -N 1
#SBATCH -t 0-02:00
#SBATCH --array=1-22
#SBATCH -p short
#SBATCH --mem=128000
#SBATCH -o /home/jz286/WES_analysis/GDReg/job_info/hostname_%j.out  
#SBATCH -e /home/jz286/WES_analysis/GDReg/job_info/hostname_%j.err 
#SBATCH --mail-type=NONE#SBATCH --mail-type=NONE


CHROM=$SLURM_ARRAY_TASK_ID
OUTPUT_PATH=/n/groups/price/martin/data_GDREG/imp_geno

# Get PLINK file from bgen
# MAF>0.1% INFO>0.6
plink2 \
    --bgen /n/groups/price/UKBiobank/download_500K/ukb_imp_chr${CHROM}_v3.bgen ref-unknown\
    --sample /n/groups/price/UKBiobank/download_500K/ukb14048_imp_chr1_v3_s487395.sample\
    --keep /n/groups/price/martin/data_GDREG/imp_geno/unrelated_337K.txt\
    --extract /n/groups/price/martin/data_GDREG/imp_geno/snp_list/snp_list_chr${CHROM}.MAF_001_INFO_06.txt\
    --rm-dup force-first\
    --maj-ref\
    --geno 0.1\
    --maf 0.001\
    --hwe 1e-50\
    --make-pgen \
    --out ${OUTPUT_PATH}/ukb_imp_chr${CHROM}_v3\
    && rm ${OUTPUT_PATH}/ukb_imp_chr${CHROM}_v3.log

# Recompute the frq info
plink2 \
    --pfile ${OUTPUT_PATH}/ukb_imp_chr${CHROM}_v3\
    --freq --out ${OUTPUT_PATH}/ukb_imp_chr${CHROM}_v3\
        && rm ${OUTPUT_PATH}/ukb_imp_chr${CHROM}_v3.log
  
# Create a small vcf file
SMALL_ID_LIST=/n/groups/price/martin/data_GDREG/imp_geno/unrelated_337K.small.txt
plink2 \
    --pfile ${OUTPUT_PATH}/ukb_imp_chr${CHROM}_v3\
    --keep $SMALL_ID_LIST --recode vcf\
    --out ${OUTPUT_PATH}/vcf/ukb_imp_chr${CHROM}_v3\
    && rm ${OUTPUT_PATH}/vcf/ukb_imp_chr${CHROM}_v3.log