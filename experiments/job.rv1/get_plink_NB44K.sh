#!/bin/bash
#SBATCH -c 4
#SBATCH -N 1
#SBATCH -t 0-12:00
#SBATCH --array=3
#SBATCH -p short
#SBATCH --mem=196000
#SBATCH -o /home/jz286/WES_analysis/LDSPEC/job_info/hostname_%j.out  
#SBATCH -e /home/jz286/WES_analysis/LDSPEC/job_info/hostname_%j.err 
#SBATCH --mail-type=NONE#SBATCH --mail-type=NONE


CHROM=$SLURM_ARRAY_TASK_ID
# CHROM=22

# Get PLINK file from bgen
OUTPUT_PATH=/n/scratch/users/j/jz286/imp_geno_NB44K
# MAF>0.1% INFO>0.6
plink2_a37 \
    --bgen /n/groups/price/UKBiobank/download_500K/ukb_imp_chr${CHROM}_v3.bgen ref-unknown\
    --sample /n/groups/price/UKBiobank/download_500K/ukb14048_imp_chr1_v3_s487395.sample\
    --keep /n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001/sample_nonBritish_White_unrelated_44K.txt\
    --extract /n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001/snp_info/snp_list_chr${CHROM}.MAF_001_INFO_06.txt\
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

# Convert to ancestral alleles (chimp): directly converting to pgen
REF_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001/snp_info/ukb_imp_chr${CHROM}_v3_aa.pvar
OUTPUT_PATH=/n/scratch/users/j/jz286/imp_geno_NB44K_chimp
plink2_a37 \
    --pfile /n/scratch/users/j/jz286/imp_geno_NB44K/ukb_imp_chr${CHROM}_v3\
    --ref-allele ${REF_FILE} 10 3\
    --make-pgen \
    --out ${OUTPUT_PATH}/ukb_imp_chr${CHROM}_v3_chimp\
    && rm ${OUTPUT_PATH}/ukb_imp_chr${CHROM}_v3_chimp.log

plink2_a37 \
    --pfile ${OUTPUT_PATH}/ukb_imp_chr${CHROM}_v3_chimp\
    --freq --out ${OUTPUT_PATH}/ukb_imp_chr${CHROM}_v3_chimp\
        && rm ${OUTPUT_PATH}/ukb_imp_chr${CHROM}_v3_chimp.log