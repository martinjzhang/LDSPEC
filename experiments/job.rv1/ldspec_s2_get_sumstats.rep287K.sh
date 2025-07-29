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

PGEN_FILE=/n/scratch/users/j/jz286/imp_geno_rep287K_chimp/ukb_imp_chr@_v3_chimp

# compute_sumstats
for i_line in {1..32}
# for i_line in 1
do
TRAIT=$( head -n $i_line "/n/groups/price/martin/LDSPEC_data/UKBB_trait_rep287K/trait_list_all_indpt.txt" | tail -1 )
PHEN_FILE=/n/groups/price/martin/LDSPEC_data/UKBB_trait_rep287K/${TRAIT}.resid.phen
PREFIX_OUT=/n/groups/price/martin/LDSPEC_data/UKBimp_rep287K_MAF001_chimp/sumstats/${TRAIT}

sbatch -p medium -t 4-00:00 -n 1 -c 4 --mem=96000 --open-mode=truncate -o $PREFIX_OUT.compute_sumstats.sbatch.log --wrap " \
python3 /home/jz286/WES_analysis/LDSPEC/CLI_simulation.py\
    --job compute_sumstats\
    --pgen_file $PGEN_FILE\
    --phen_file $PHEN_FILE\
    --random_seed 0\
    --prefix_out $PREFIX_OUT"
done
