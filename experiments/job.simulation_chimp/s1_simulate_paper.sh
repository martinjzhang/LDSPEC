#!/bin/bash

PGEN_FILE=/n/scratch3/users/j/jz286/imp_geno_chimp/ukb_imp_chr1_v3_chimp
ANNOT_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.simulation_chimp/reg_annot_file.sim.txt

# 128gb
SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.100123/causal_h2g50_p20
SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.100123/causal_h2g50_p10
SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.100123/causal_h2g20_p20
SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.100123/causal_pos_h2g50_p20
SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.100123/causal_neg_h2g50_p20
# # 32gb
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.100123/null_h2g50_p20
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.100123/null_h2g50_p10
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.100123/null_h2g20_p20


for i_rep in {0..49}
# for i_rep in 0
do
CONFIG_FILE=${SIMU_PATH}/config
PREFIX_OUT=${SIMU_PATH}/rep${i_rep}
    
# # Independent sparsity
# sbatch -p medium -t 0-24:00 -n 1 -c 1 --mem=32000 --open-mode=truncate -o $PREFIX_OUT.sbatch.log --wrap " \
# python3 /home/jz286/WES_analysis/LDSPEC/CLI_simulation.py\
#     --job simulate\
#     --pgen_file $PGEN_FILE\
#     --config_file $CONFIG_FILE\
#     --annot_file $ANNOT_FILE\
#     --random_seed $i_rep\
#     --prefix_out $PREFIX_OUT"
    
# Block-wise sparsity
sbatch -p medium -t 0-24:00 -n 1 -c 1 --mem=128000 --open-mode=truncate -o $PREFIX_OUT.sbatch.log --wrap " \
python3 /home/jz286/WES_analysis/LDSPEC/CLI_simulation.py\
    --job simulate\
    --pgen_file $PGEN_FILE\
    --config_file $CONFIG_FILE\
    --annot_file $ANNOT_FILE\
    --random_seed $i_rep\
    --flag_bw_sparse True\
    --prefix_out $PREFIX_OUT"
done