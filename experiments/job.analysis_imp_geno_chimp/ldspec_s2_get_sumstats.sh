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

PGEN_FILE=/n/scratch3/users/j/jz286/imp_geno_chimp/ukb_imp_chr@_v3_chimp

# compute_sumstats
for i_line in {2..118}
# for i_line in 1
do
TRAIT=$( head -n $i_line "/n/groups/price/martin/LDSPEC_data/UKBB_trait/trait_list_all.txt" | tail -1 )
PHEN_FILE=/n/groups/price/martin/LDSPEC_data/UKBB_trait/${TRAIT}.resid.phen
PREFIX_OUT=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/sumstats/${TRAIT}

sbatch -p medium -t 4-00:00 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $PREFIX_OUT.compute_sumstats.sbatch.log --wrap " \
python3 /home/jz286/WES_analysis/LDSPEC/CLI_simulation.py\
    --job compute_sumstats\
    --pgen_file $PGEN_FILE\
    --phen_file $PHEN_FILE\
    --random_seed 0\
    --prefix_out $PREFIX_OUT"
done
