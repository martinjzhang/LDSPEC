#!/bin/bash

PGEN_FILE=/n/scratch/users/j/jz286/imp_geno_chimp/ukb_imp_chr1_v3_chimp
ANNOT_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.simulation_chimp/reg_annot_file.sim.rv1.txt

# 128gb
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.011425/causal_h2g50_p20
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.011425/causal_h2g50_p1
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.011425/causal_h2g50_p100
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.011425/causal_highld_h2g50_p20
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.011425/causal_mafdepend_h2g50_p20
SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.011425/causal_simple_h2g50_p20
# 32gb
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.011425/null_h2g50_p20
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.011425/null_h2g50_p1
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.011425/null_h2g50_p100


for i_rep in {1..49}
# for i_rep in 0
do
CONFIG_FILE=${SIMU_PATH}/config
PREFIX_OUT=${SIMU_PATH}/rep${i_rep}
    
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



# # simulation using rep287K
# PGEN_FILE=/n/scratch/users/j/jz286/imp_geno_rep287K_chimp/ukb_imp_chr1_v3_chimp
# ANNOT_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.simulation_chimp/reg_annot_file.sim.rep287K.rv1.txt

# # 128gb
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.011425/causal_h2g50_p20.rep287K
# # 32gb
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.011425/null_h2g50_p20.rep287K


# # for i_rep in {1..49}
# for i_rep in 0
# do
# CONFIG_FILE=${SIMU_PATH}/config
# PREFIX_OUT=${SIMU_PATH}/rep${i_rep}
    
# # Block-wise sparsity
# sbatch -p medium -t 0-24:00 -n 1 -c 1 --mem=32000 --open-mode=truncate -o $PREFIX_OUT.sbatch.log --wrap " \
# python3 /home/jz286/WES_analysis/LDSPEC/CLI_simulation.py\
#     --job simulate\
#     --pgen_file $PGEN_FILE\
#     --config_file $CONFIG_FILE\
#     --annot_file $ANNOT_FILE\
#     --random_seed $i_rep\
#     --flag_bw_sparse True\
#     --prefix_out $PREFIX_OUT"
# done