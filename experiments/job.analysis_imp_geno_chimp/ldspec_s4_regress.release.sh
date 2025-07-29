#!/bin/bash

LDSPEC_DATA_PATH=/n/groups/price/martin/LDSPEC_data/LDSPEC_release_100723
PGEN_FILE=/n/scratch3/users/j/jz286/imp_geno_chimp/ukb_imp_chr@_v3_chimp
AVGR_FILE=${LDSPEC_DATA_PATH}/main.avgr

ANALYSIS=prox_gene_fct_all_ld
ANNOT_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.analysis_imp_geno_chimp/reg_annot_file/reg_annot_file.baseline_sp_release.txt
SCORE_FILE=${LDSPEC_DATA_PATH}/score_single165/baseline.c@_score.tsv.gz,${LDSPEC_DATA_PATH}/score_pair136/baseline_sp.c@_score.tsv.gz

for i_line in 1
do
TRAIT=$( head -n $i_line ${LDSPEC_DATA_PATH}/trait_list_all.txt | tail -1 ) # 29 indpt analyzed in the paper

SUMSTATS_FILE=${LDSPEC_DATA_PATH}/sumstats/${TRAIT}.nomhc.sumstats.gz
PATH_OUT=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_res_092223.${ANALYSIS}.release
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