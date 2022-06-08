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

PGEN_FILE=/n/groups/price/martin/WES_analysis/toy_1K/chr@_v1.SPB.hg19.toy_1K

# Full LD matrix
for CHR in {1..10}
# for CHR in 1
do
PREFIX_OUT=/n/groups/price/martin/WES_analysis/toy_1K/results/top_1K_chr${CHR}
python3 /home/jz286/WES_analysis/GDReg/run_gdreg.py\
    --job compute_ld\
    --pgen_file $PGEN_FILE\
    --snp_range 'chr='${CHR}'|chr_ref='${CHR}\
    --random_seed 0\
    --memory 2048\
    --flag_full_ld True\
    --prefix_out $PREFIX_OUT
done

# Reduce LD matrix (defaut)
for CHR in {1..10}
# for CHR in 1
do
PREFIX_OUT=/n/groups/price/martin/WES_analysis/toy_1K/results/top_1K_chr${CHR}
python3 /home/jz286/WES_analysis/GDReg/run_gdreg.py\
    --job compute_ld\
    --pgen_file $PGEN_FILE\
    --snp_range 'chr='${CHR}'|chr_ref='${CHR}\
    --random_seed 0\
    --memory 2048\
    --prefix_out $PREFIX_OUT
done


# for CHR in {1..10}
# do
# for CHR_REF in {1..10}
# do
# PREFIX_OUT=/n/groups/price/martin/WES_analysis/toy_1K/results/full_ld/top_1K_chr${CHR}_chr${CHR_REF}
# python3 /home/jz286/WES_analysis/GDReg/run_gdreg.py\
#     --job compute_ld\
#     --pgen_file $PGEN_FILE\
#     --snp_range 'chr='${CHR}'|chr_ref='${CHR_REF}\
#     --random_seed 0\
#     --memory 256\
#     --prefix_out $PREFIX_OUT
# done
# done
