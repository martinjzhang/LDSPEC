#!/bin/bash
#SBATCH -c 8
#SBATCH -N 1
#SBATCH -t 0-01:00
#SBATCH --array=1-7
#SBATCH -p short
#SBATCH --mem=64000
#SBATCH -o /home/jz286/WES_analysis/GDReg/experiments/job_info/hostname_%j.out  
#SBATCH -e /home/jz286/WES_analysis/GDReg/experiments/job_info/hostname_%j.err 
#SBATCH --mail-type=NONE#SBATCH --mail-type=NONE

PGEN_FILE=/n/scratch3/users/j/jz286/imp_geno_chimp/ukb_imp_chr@_v3_chimp

# # get_snp_block
# PREFIX_OUT=/n/scratch3/users/j/jz286/imp_geno.gdreg_ld/ukb_imp_v3
# python3 /home/jz286/WES_analysis/GDReg/run_gdreg.py\
#     --job get_snp_block\
#     --pgen_file $PGEN_FILE\
#     --prefix_out $PREFIX_OUT
    
    
# compute_ld
# chr1: 117 files
for i_line in {1..2}
# for i_line in 1
do 
PREFIX_OUT=/n/scratch3/users/j/jz286/imp_geno_chimp.gdreg_ld_win1e7/ukb_imp_v3_chimp
SNP_RANGE=$( head -n $i_line "/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/ukb_imp_v3.snp_range.uf.txt" | tail -1 )

# python3 /home/jz286/WES_analysis/GDReg/run_gdreg.py\
#     --job compute_ld\
#     --pgen_file $PGEN_FILE\
#     --snp_range $SNP_RANGE\
#     --win_size 1e6\
#     --prefix_out $PREFIX_OUT

sbatch -p medium -t 0-36:00 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $PREFIX_OUT.${SNP_RANGE}.compute_ld.sbatch.log --wrap " \
python3 /home/jz286/WES_analysis/GDReg/run_gdreg.py\
    --job compute_ld\
    --pgen_file $PGEN_FILE\
    --snp_range $SNP_RANGE\
    --win_size 1e7\
    --prefix_out $PREFIX_OUT"
done