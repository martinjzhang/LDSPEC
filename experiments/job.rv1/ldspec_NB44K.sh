#!/bin/bash
#SBATCH -c 8
#SBATCH -N 1
#SBATCH -t 0-01:00
#SBATCH --array=1
#SBATCH -p short
#SBATCH --mem=64000
#SBATCH -o /home/jz286/WES_analysis/GDReg/experiments/job_info/hostname_%j.out  
#SBATCH -e /home/jz286/WES_analysis/GDReg/experiments/job_info/hostname_%j.err 
#SBATCH --mail-type=NONE#SBATCH --mail-type=NONE

PGEN_FILE=/n/scratch/users/j/jz286/imp_geno_NB44K_chimp/ukb_imp_chr@_v3_chimp

# # get_snp_block
# PREFIX_OUT=/n/groups/price/martin/LDSPEC_data/UKBimp_NK44K_MAF001_chimp/ukb_imp_v3
# python3 /home/jz286/WES_analysis/LDSPEC/CLI_ldspec.py\
#     --job get_snp_block\
#     --pgen_file $PGEN_FILE\
#     --prefix_out $PREFIX_OUT
    
    
# # compute_ld (chr1: 121 files)
# # for i_line in {3..1439}
# for i_line in {1..111}
# do 
# PREFIX_OUT=/n/scratch/users/j/jz286/imp_geno_NB44K_chimp.ldspec_ld_1e7/ukb_imp_v3_chimp
# SNP_RANGE=$( head -n $i_line "/n/groups/price/martin/LDSPEC_data/UKBimp_NB44K_MAF001_chimp/ukb_imp_v3.snp_range.txt" | tail -1 )
# SNP_RANGE=$( head -n $i_line "/n/groups/price/martin/LDSPEC_data/UKBimp_NB44K_MAF001_chimp/ukb_imp_v3.snp_range.uf.txt" | tail -1 )

# # Initial: sbatch -p short -t 00-6:00 -n 1 -c 1 --mem=16000
# # Second: sbatch -p short -t 00-12:00 -n 1 -c 1 --mem=24000
# sbatch -p short -t 00-6:00 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $PREFIX_OUT.${SNP_RANGE}.compute_ld.sbatch.log --wrap " \
# python3 /home/jz286/WES_analysis/LDSPEC/CLI_ldspec.py\
#     --job compute_ld\
#     --pgen_file $PGEN_FILE\
#     --snp_range $SNP_RANGE\
#     --win_size 1e7\
#     --prefix_out $PREFIX_OUT"
# done


# # compute_score
# for i_line in {1..111}
# # for i_line in {1..2}
# do
# SNP_RANGE=$( head -n $i_line '/n/groups/price/martin/LDSPEC_data/UKBimp_NB44K_MAF001_chimp/ukb_imp_v3.snp_range.txt' | tail -1 )
# SNP_RANGE=$( head -n $i_line '/n/groups/price/martin/LDSPEC_data/UKBimp_NB44K_MAF001_chimp/ukb_imp_v3.snp_range.uf.txt' | tail -1 )
# LD_FILE=/n/scratch/users/j/jz286/imp_geno_NB44K_chimp.ldspec_ld_1e7/ukb_imp_v3_chimp.${SNP_RANGE}_ld.npz

# ANNOT_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.analysis_imp_geno_chimp/reg_annot_file/reg_annot_file.NB44K.prox_gene_fct_all_ld.txt
# PREFIX_OUT=/n/scratch/users/j/jz286/imp_geno_NB44K_chimp.score/ukb_imp_v3
    
# sbatch -p short -t 0-2:00 -n 1 -c 1 --mem=24000 --open-mode=truncate -o ${PREFIX_OUT}.${SNP_RANGE}.compute_score.sbatch.log --wrap " \
# python3 /home/jz286/WES_analysis/LDSPEC/CLI_ldspec.py\
#     --job compute_score\
#     --pgen_file $PGEN_FILE\
#     --ld_file $LD_FILE\
#     --annot_file $ANNOT_FILE\
#     --win_size 1e7\
#     --prefix_out $PREFIX_OUT"
# done


# # combine_score
# SCORE_FILE=/n/scratch/users/j/jz286/imp_geno_NB44K_chimp.score/ukb_imp_v3.@_score.tsv.gz
# SNP_RANGE_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_NB44K_MAF001_chimp/ukb_imp_v3.snp_range.txt
# PREFIX_OUT=/n/groups/price/martin/LDSPEC_data/UKBimp_NB44K_MAF001_chimp/baseline_sp
# sbatch -p short -t 0-12:00 -n 1 -c 1 --mem=96000 --open-mode=truncate -o ${PREFIX_OUT}.${SNP_RANGE}.compute_score.sbatch.log --wrap " \
# python3 /home/jz286/WES_analysis/LDSPEC/CLI_ldspec.py\
#     --job combine_score\
#     --score_file $SCORE_FILE\
#     --snp_range_file $SNP_RANGE_FILE\
#     --prefix_out $PREFIX_OUT"
    


# # compute_avgr
# for i_line in {2..137}
# do 
# ANNOT_FILE=$( head -n $i_line "/home/jz286/WES_analysis/LDSPEC/experiments/job.analysis_imp_geno_chimp/reg_annot_file/reg_annot_file.NB44K.prox_gene_fct_all_ld.txt" | tail -1 )
# PREFIX_OUT=$ANNOT_FILE
# LD_FILE=/n/scratch/users/j/jz286/imp_geno_NB44K_chimp.ldspec_ld_1e7/ukb_imp_v3_chimp.@_ld.npz
    
# sbatch -p short -t 0-3:00 -n 1 -c 1 --mem=32000 --open-mode=truncate -o ${PREFIX_OUT}.compute_avgr.sbatch.log --wrap " \
# python3 /home/jz286/WES_analysis/LDSPEC/CLI_ldspec.py\
#     --job compute_avgr\
#     --pgen_file $PGEN_FILE\
#     --ld_file $LD_FILE\
#     --annot_file $ANNOT_FILE\
#     --prefix_out $PREFIX_OUT"
# done


# regress
ANALYSIS=prox_gene_fct_all_ld_NB44K
ANNOT_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.analysis_imp_geno_chimp/reg_annot_file/reg_annot_file.NB44K.prox_gene_fct_all_ld.txt
SCORE_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_NB44K_MAF001_chimp/baseline_sp.c@_score.tsv.gz
AVGR_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_NB44K_MAF001_chimp/main.avgr

# for i_line in 1
for i_line in {2..29}
do
TRAIT=$( head -n $i_line /n/groups/price/martin/LDSPEC_data/UKBB_trait/trait_list_paper_indpt.txt | tail -1 ) # 29 indpt analyzed in the paper

SUMSTATS_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_rep287K_MAF001_chimp/sumstats/${TRAIT}.nomhc.sumstats.gz
PATH_OUT=/n/groups/price/martin/LDSPEC_data/UKBimp_rep287K_MAF001_chimp/ldspec_res_040925.${ANALYSIS}
PREFIX_OUT=${PATH_OUT}/${TRAIT}
[ -d $PATH_OUT ] || mkdir $PATH_OUT
    
sbatch -p medium -t 0-24:00 -n 1 -c 1 --mem=128000 --open-mode=truncate -o $PREFIX_OUT.sbatch.log --wrap " \
python3 /home/jz286/WES_analysis/LDSPEC/CLI_ldspec.py\
    --job regress\
    --pgen_file $PGEN_FILE\
    --annot_file $ANNOT_FILE\
    --score_file $SCORE_FILE\
    --sumstats_file $SUMSTATS_FILE\
    --avgr_file $AVGR_FILE\
    --prefix_out $PREFIX_OUT"
done