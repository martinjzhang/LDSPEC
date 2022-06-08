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

# # compute_phen
# for i_rep in {0..19}
# do

# EFF_FILE=/n/groups/price/martin/WES_analysis/toy_1K/sanity_nd_rep${i_rep}.eff.gz
# PREFIX_OUT=/n/groups/price/martin/WES_analysis/toy_1K/sanity_nd_rep${i_rep}

# sbatch -p short -t 0-00:15 -n 1 -c 1 --mem=16000 --open-mode=truncate -e $PREFIX_OUT.sbatch.err -o $PREFIX_OUT.sbatch.log --wrap " \
# python3 /home/jz286/WES_analysis/GDReg/run_simulation.py\
#     --job compute_phen\
#     --pgen_file $PGEN_FILE\
#     --eff_file $EFF_FILE\
#     --random_seed 0\
#     --prefix_out $PREFIX_OUT"
# done

# compute_sumstats
for i_rep in {0..19}
do
PHEN_FILE=/n/groups/price/martin/WES_analysis/toy_1K/sanity_nd_rep${i_rep}.phen
PREFIX_OUT=/n/groups/price/martin/WES_analysis/toy_1K/sanity_nd_rep${i_rep}

sbatch -p short -t 0-00:15 -n 1 -c 1 --mem=16000 --open-mode=truncate -e $PREFIX_OUT.sbatch.err -o $PREFIX_OUT.sbatch.log --wrap " \
python3 /home/jz286/WES_analysis/GDReg/run_simulation.py\
    --job compute_sumstats\
    --pgen_file $PGEN_FILE\
    --phen_file $PHEN_FILE\
    --random_seed 0\
    --prefix_out $PREFIX_OUT"
done
