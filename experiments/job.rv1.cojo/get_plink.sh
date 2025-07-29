#!/bin/bash
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-2:00
#SBATCH --array=1-22
#SBATCH -p short
#SBATCH --mem=128000
#SBATCH -o /home/jz286/WES_analysis/LDSPEC/job_info/hostname_%j.out  
#SBATCH -e /home/jz286/WES_analysis/LDSPEC/job_info/hostname_%j.err 
#SBATCH --mail-type=NONE#SBATCH --mail-type=NONE


CHROM=$SLURM_ARRAY_TASK_ID
# CHROM=22

OUTPUT_PATH=/n/scratch/users/j/jz286/imp_geno_chimp_bed_MAF01
plink2_a37 \
    --pfile /n/scratch/users/j/jz286/imp_geno_chimp/ukb_imp_chr${CHROM}_v3_chimp\
    --maf 0.01\
    --make-bed\
    --out ${OUTPUT_PATH}/ukb_imp_chr${CHROM}_v3_chimp\
        && rm ${OUTPUT_PATH}/ukb_imp_chr${CHROM}_v3_chimp.log
        
plink2 \
    --bfile ${OUTPUT_PATH}/ukb_imp_chr${CHROM}_v3_chimp\
    --freq --out ${OUTPUT_PATH}/ukb_imp_chr${CHROM}_v3_chimp\
        && rm ${OUTPUT_PATH}/ukb_imp_chr${CHROM}_v3_chimp.log