#!/bin/bash

PGEN_FILE=/n/scratch3/users/j/jz286/imp_geno_chimp/ukb_imp_chr@_v3_chimp
AVGR_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot/main.avgr

# ANALYSIS=bsl
# ANNOT_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001/baseline_annot/baseline_165annots_chr@.annot.gz
# SCORE_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/baseline.c@_score.tsv.gz

# ANALYSIS=prox_gene_fct_all_ld
# ANNOT_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.analysis_imp_geno_chimp/reg_annot_file/reg_annot_file.prox_gene_fct_all_ld.txt
# SCORE_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/baseline.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/pannot_prox_ld.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/pannot_gene_ld.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/pannot_baseline_0_100_ld.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/pannot_baseline_0_1000_ld.c@_score.tsv.gz

ANALYSIS=prox_all_ld
ANNOT_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.analysis_imp_geno_chimp/reg_annot_file/reg_annot_file.prox_all_ld.txt
SCORE_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/baseline.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/pannot_prox_ld.c@_score.tsv.gz

# ANALYSIS=gene_all_ld
# ANNOT_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.analysis_imp_geno_chimp/reg_annot_file/reg_annot_file.gene_all_ld.txt
# SCORE_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/baseline.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/pannot_gene_ld.c@_score.tsv.gz

# ANALYSIS=fct_all_ld
# ANNOT_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.analysis_imp_geno_chimp/reg_annot_file/reg_annot_file.fct_all_ld.txt
# SCORE_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/baseline.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/pannot_baseline_0_100_ld.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/pannot_baseline_0_1000_ld.c@_score.tsv.gz

for i_line in {2..29}
# for i_line in 1
do
# TRAIT=$( head -n $i_line "/n/groups/price/martin/LDSPEC_data/UKBB_trait/trait_list_all.txt" | tail -1 )
TRAIT=$( head -n $i_line "/n/groups/price/martin/LDSPEC_data/UKBB_trait/trait_list_paper_indpt.txt" | tail -1 ) # 29 indpt analyzed in the paper
# TRAIT=$( head -n $i_line "/n/groups/price/martin/LDSPEC_data/UKBB_trait/trait_list_paper_other.txt" | tail -1 ) # 29 indpt analyzed in the paper
SUMSTATS_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/sumstats/${TRAIT}.nomhc.sumstats.gz
PATH_OUT=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_res_092223.${ANALYSIS}
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