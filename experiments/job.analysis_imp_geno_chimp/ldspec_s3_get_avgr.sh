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
for i_line in {1..76}
# for i_line in 1
do 
ANNOT_FILE=$( head -n $i_line "/home/jz286/WES_analysis/LDSPEC/experiments/job.curate_data/pannot_list/pannot_list_basic.txt" | tail -1 )
ANNOT_FILE=$( head -n $i_line "/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot/missing_avgr.txt" | tail -1 )
ANNOT_FILE=$( head -n $i_line "/home/jz286/WES_analysis/LDSPEC/experiments/job.curate_data/pannot_list/pannot_list_gene.txt" | tail -1 )
ANNOT_FILE=$( head -n $i_line "/home/jz286/WES_analysis/LDSPEC/experiments/job.curate_data/pannot_list/pannot_list_baseline_0_100_ld.txt" | tail -1 )
ANNOT_FILE=$( head -n $i_line "/home/jz286/WES_analysis/LDSPEC/experiments/job.curate_data/pannot_list/pannot_list_baseline_0_1000_ld.txt" | tail -1 )
PREFIX_OUT=$ANNOT_FILE
LD_FILE=/n/scratch3/users/j/jz286/imp_geno_chimp.gdreg_ld_1e7/ukb_imp_v3_chimp.@_ld.npz
    
sbatch -p short -t 0-8:00 -n 1 -c 1 --mem=32000 --open-mode=truncate -o ${PREFIX_OUT}.compute_avgr.sbatch.log --wrap " \
python3 /home/jz286/WES_analysis/LDSPEC/CLI_ldspec.py\
    --job compute_avgr\
    --pgen_file $PGEN_FILE\
    --ld_file $LD_FILE\
    --annot_file $ANNOT_FILE\
    --prefix_out $PREFIX_OUT"
done