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

# # Get a pgen file
# plink2 \
#     --bfile /n/groups/price/martin/WES_analysis/WES_34K_maf1en4/chr20_v1.SPB.hg19 \
#     --geno 0.1 \
#     --mac 1 \
#     --hwe 1e-50\
#     --maj-ref \
#     --make-pgen \
#     --out /n/groups/price/martin/WES_analysis/toy_1K/results_wes/chr20_v1.SPB.hg19

# # Recompute the frq info
# plink2 \
#     --pfile /n/groups/price/martin/WES_analysis/toy_1K/results_wes/chr20_v1.SPB.hg19\
#     --freq \
#     --out /n/groups/price/martin/WES_analysis/toy_1K/results_wes/chr20_v1.SPB.hg19


PGEN_FILE=/n/groups/price/martin/WES_analysis/toy_1K/results_wes/chr@_v1.SPB.hg19
PREFIX_OUT=/n/groups/price/martin/WES_analysis/toy_1K/results_wes/WES_34K_chr20
python3 /home/jz286/WES_analysis/GDReg/run_gdreg.py\
    --job compute_ld\
    --pgen_file $PGEN_FILE\
    --snp_range 'chr=20|chr_ref=20'\
    --random_seed 0\
    --memory 8024\
    --prefix_out $PREFIX_OUT