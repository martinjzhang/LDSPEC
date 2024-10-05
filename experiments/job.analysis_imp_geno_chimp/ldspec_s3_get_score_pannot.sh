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

PGEN_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/imp_geno_chimp/ukb_imp_chr@_v3_chimp
PGEN_FILE=/n/scratch3/users/j/jz286/imp_geno_chimp/ukb_imp_chr@_v3_chimp

# for i_line in {2..1492}
for i_line in 1492
do
SNP_RANGE=$( head -n $i_line '/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ukb_imp_v3.snp_range.txt' | tail -1 )
LD_FILE=/n/scratch3/users/j/jz286/imp_geno_chimp.ldspec_ld_1e7/ukb_imp_v3_chimp.${SNP_RANGE}_ld.npz

ANNOT_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001/baseline_annot/baseline_165annots_chr@.annot.gz
PREFIX_OUT=/n/scratch3/users/j/jz286/imp_geno_chimp.score_annot/ukb_imp_v3
    
sbatch -p short -t 0-02:00 -n 1 -c 1 --mem=24000 --open-mode=truncate -o ${PREFIX_OUT}.${SNP_RANGE}.compute_score.sbatch.log --wrap " \
python3 /home/jz286/WES_analysis/LDSPEC/CLI_ldspec.py\
    --job compute_score\
    --pgen_file $PGEN_FILE\
    --ld_file $LD_FILE\
    --annot_file $ANNOT_FILE\
    --win_size 1e7\
    --prefix_out $PREFIX_OUT"

# ANNOT_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.curate_data/pannot_list/pannot_list_gene.txt
# PREFIX_OUT=/n/scratch3/users/j/jz286/imp_geno_chimp.score_pannot_gene/ukb_imp_v3
# sbatch -p short -t 0-02:00 -n 1 -c 1 --mem=24000 --open-mode=truncate -o ${PREFIX_OUT}.${SNP_RANGE}.compute_score.sbatch.log --wrap " \
# python3 /home/jz286/WES_analysis/LDSPEC/CLI_ldspec.py\
#     --job compute_score\
#     --pgen_file $PGEN_FILE\
#     --ld_file $LD_FILE\
#     --annot_file $ANNOT_FILE\
#     --win_size 1e7\
#     --prefix_out $PREFIX_OUT"

# ANNOT_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.curate_data/pannot_list/pannot_list_baseline_0_100_ld.txt
# PREFIX_OUT=/n/scratch3/users/j/jz286/imp_geno_chimp.score_pannot_baseline_0_100_ld/ukb_imp_v3
# sbatch -p short -t 0-01:00 -n 1 -c 1 --mem=24000 --open-mode=truncate -o ${PREFIX_OUT}.${SNP_RANGE}.compute_score.sbatch.log --wrap " \
# python3 /home/jz286/WES_analysis/LDSPEC/CLI_ldspec.py\
#     --job compute_score\
#     --pgen_file $PGEN_FILE\
#     --ld_file $LD_FILE\
#     --annot_file $ANNOT_FILE\
#     --win_size 1e7\
#     --prefix_out $PREFIX_OUT"
    
# ANNOT_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.curate_data/pannot_list/pannot_list_baseline_0_1000_ld.txt
# PREFIX_OUT=/n/scratch3/users/j/jz286/imp_geno_chimp.score_pannot_baseline_0_1000_ld/ukb_imp_v3
# sbatch -p short -t 0-01:00 -n 1 -c 1 --mem=24000 --open-mode=truncate -o ${PREFIX_OUT}.${SNP_RANGE}.compute_score.sbatch.log --wrap " \
# python3 /home/jz286/WES_analysis/LDSPEC/CLI_ldspec.py\
#     --job compute_score\
#     --pgen_file $PGEN_FILE\
#     --ld_file $LD_FILE\
#     --annot_file $ANNOT_FILE\
#     --win_size 1e7\
#     --prefix_out $PREFIX_OUT"
done