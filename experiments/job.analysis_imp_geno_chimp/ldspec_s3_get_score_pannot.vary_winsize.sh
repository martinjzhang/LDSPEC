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

# chr1: 117 files
# for WIN_SIZE in {1e6,3e6,5e6,1e7}
for WIN_SIZE in {1e6,3e6,5e6}
do
for i_line in {2..117}
# for i_line in 1
do
SNP_RANGE=$( head -n $i_line "/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ukb_imp_v3.snp_range.txt" | tail -1 )
LD_FILE=/n/scratch3/users/j/jz286/imp_geno_chimp.gdreg_ld_1e7/ukb_imp_v3_chimp.${SNP_RANGE}_ld.npz

ANNOT_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.curate_data/pannot_list/pannot_list_basic.txt
PATH_OUT=/n/scratch3/users/j/jz286/imp_geno_chimp.score_pannot_basic.${WIN_SIZE}

# ANNOT_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.curate_data/pannot_list/pannot_list_gene.txt
# PATH_OUT=/n/scratch3/users/j/jz286/imp_geno_chimp.score_pannot_gene.${WIN_SIZE}

# ANNOT_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.curate_data/pannot_list/pannot_list_baseline_0_100_ld.txt
# PATH_OUT=/n/scratch3/users/j/jz286/imp_geno_chimp.score_pannot_baseline_0_100_ld.${WIN_SIZE}

# ANNOT_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.curate_data/pannot_list/pannot_list_baseline_0_1000_ld.txt
# PATH_OUT=/n/scratch3/users/j/jz286/imp_geno_chimp.score_pannot_baseline_0_1000_ld.${WIN_SIZE}

[ -d $PATH_OUT ] || mkdir $PATH_OUT
PREFIX_OUT=${PATH_OUT}/ukb_imp_v3
sbatch -p short -t 0-01:00 -n 1 -c 1 --mem=24000 --open-mode=truncate -o ${PREFIX_OUT}.${SNP_RANGE}.compute_score.sbatch.log --wrap " \
python3 /home/jz286/WES_analysis/LDSPEC/CLI_ldspec.py\
    --job compute_score\
    --pgen_file $PGEN_FILE\
    --ld_file $LD_FILE\
    --annot_file $ANNOT_FILE\
    --win_size $WIN_SIZE\
    --prefix_out $PREFIX_OUT"
done
done