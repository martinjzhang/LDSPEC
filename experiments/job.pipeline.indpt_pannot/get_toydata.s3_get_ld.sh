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

PGEN_FILE=/n/groups/price/martin/WES_analysis/toy_10K/chr@_v1.SPB.hg19.toy_10K

# # get_snp_block
# PREFIX_OUT=/n/groups/price/martin/WES_analysis/toy_10K/res_gdreg/toy_10K
# python3 /home/jz286/WES_analysis/GDReg/run_gdreg.py\
#     --job get_snp_block\
#     --pgen_file $PGEN_FILE\
#     --prefix_out $PREFIX_OUT
    
    
# # compute_ld
# for i_line in {1..20}
# do 
# PREFIX_OUT=/n/groups/price/martin/WES_analysis/toy_10K/res_gdreg/toy_10K
# SNP_RANGE=$( head -n $i_line "/n/groups/price/martin/WES_analysis/toy_10K/res_gdreg/toy_10K.snp_range.txt" | tail -1 )
# sbatch -p short -t 0-00:15 -n 1 -c 1 --mem=8000 --open-mode=truncate -o $PREFIX_OUT.compute_ld.sbatch.log --wrap " \
# python3 /home/jz286/WES_analysis/GDReg/run_gdreg.py\
#     --job compute_ld\
#     --pgen_file $PGEN_FILE\
#     --snp_range $SNP_RANGE\
#     --random_seed 0\
#     --memory 2048\
#     --prefix_out $PREFIX_OUT"
# done

# compute_ld : full_ld
for CHR in {1..2}
do 
for CHR_REF in {1..10}
do 
PREFIX_OUT=/n/groups/price/martin/WES_analysis/toy_10K/gdreg_file_fullld/toy_10K
SNP_RANGE=c${CHR}_r${CHR_REF}
sbatch -p short -t 0-00:05 -n 1 -c 1 --mem=4000 --open-mode=truncate --wrap " \
python3 /home/jz286/WES_analysis/GDReg/run_gdreg.py\
    --job compute_ld\
    --pgen_file $PGEN_FILE\
    --snp_range $SNP_RANGE\
    --flag_full_ld True\
    --random_seed 0\
    --memory 2048\
    --prefix_out $PREFIX_OUT"
done
done