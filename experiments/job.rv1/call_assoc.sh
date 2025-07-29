#!/bin/bash
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-02:00
#SBATCH --array=1
#SBATCH -p short
#SBATCH --mem=16000
#SBATCH -o /home/jz286/WES_analysis/GDReg/experiments/job_info/hostname_%j.out  
#SBATCH -e /home/jz286/WES_analysis/GDReg/experiments/job_info/hostname_%j.err 
#SBATCH --mail-type=NONE#SBATCH --mail-type=NONE
    
    
# for i_line in {1..32}
for i_line in 8
do 
TRAIT=$( head -n $i_line "/n/groups/price/martin/LDSPEC_data/UKBB_trait_ss50K/trait_list_all_indpt.txt" | tail -1 )
PREFIX_OUT=/n/groups/price/martin/LDSPEC_data/rv1_assoc/jnt_test_120724/UKBimp_ss50K_MAF001_chimp
echo ${TRAIT}

# python3 /home/jz286/WES_analysis/LDSPEC/experiments/job.rv1_assoc/assoc.py\
#     --trait ${TRAIT}\
#     --prefix_out $PREFIX_OUT

sbatch -p medium -t 0-48:00 -n 1 -c 1 --mem=64000 --open-mode=truncate -o $PREFIX_OUT.${TRAIT}.log --wrap " \
python3 /home/jz286/WES_analysis/LDSPEC/experiments/job.rv1_assoc/assoc.py\
    --trait ${TRAIT}\
    --prefix_out $PREFIX_OUT"
done