# #!/bin/bash

# original 334K simulations (simulation.100123)
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.100123/null_h2g50_p20
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.100123/causal_h2g50_p20

# original 334K simulations
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.011425/causal_h2g50_p20
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.011425/causal_h2g50_p1
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.011425/causal_h2g50_p100
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.011425/causal_highld_h2g50_p20
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.011425/causal_mafdepend_h2g50_p20
SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.011425/causal_simple_h2g50_p20
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.011425/null_h2g50_p20
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.011425/null_h2g50_p1
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.011425/null_h2g50_p100

# rep287K simulations
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.011425/causal_h2g50_p20.rep287K
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.011425/null_h2g50_p20.rep287K



ANALYSIS=prox_gene_fct_all_ld
PGEN_FILE=/n/scratch/users/j/jz286/imp_geno_chimp/ukb_imp_chr1_v3_chimp
ANNOT_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.analysis_imp_geno_chimp/reg_annot_file/reg_annot_file.prox_gene_fct_all_ld.txt
AVGR_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot/main.avgr
SCORE_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/baseline.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/pannot_prox_ld.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/pannot_gene_ld.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/pannot_baseline_0_100_ld.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/pannot_baseline_0_1000_ld.c@_score.tsv.gz


# ANALYSIS=simu_causal_pannot
# PGEN_FILE=/n/scratch/users/j/jz286/imp_geno_chimp/ukb_imp_chr1_v3_chimp
# ANNOT_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.analysis_imp_geno_chimp/reg_annot_file/reg_annot_file.simu_causal_pannot.txt
# AVGR_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot/main.avgr
# SCORE_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/baseline.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/pannot_prox_ld.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/pannot_baseline_0_1000_ld.c@_score.tsv.gz


# ANALYSIS=prox_gene_fct_all_ld_ss50K
# PGEN_FILE=/n/scratch/users/j/jz286/imp_geno_ss50K_chimp/ukb_imp_chr1_v3_chimp
# ANNOT_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.analysis_imp_geno_chimp/reg_annot_file/reg_annot_file.ss50K.prox_gene_fct_all_ld.txt
# SCORE_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_ss50K_MAF001_chimp/baseline_sp.c1_score.tsv.gz
# AVGR_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_ss50K_MAF001_chimp/main.avgr


# ANALYSIS=prox_gene_fct_all_ld_NB44K
# PGEN_FILE=/n/scratch/users/j/jz286/imp_geno_NB44K_chimp/ukb_imp_chr1_v3_chimp
# ANNOT_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.analysis_imp_geno_chimp/reg_annot_file/reg_annot_file.NB44K.prox_gene_fct_all_ld.txt
# SCORE_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_NB44K_MAF001_chimp/baseline_sp.c1_score.tsv.gz
# AVGR_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_NB44K_MAF001_chimp/main.avgr

OUT_PATH=${SIMU_PATH}.${ANALYSIS}

[ -d $OUT_PATH ] || mkdir $OUT_PATH

for i_rep in {0..49}
# for i_rep in 0
do
SUMSTATS_FILE=${SIMU_PATH}/rep${i_rep}.sumstats.gz
PREFIX_OUT=${OUT_PATH}/rep${i_rep}

sbatch -p short -t 0-3:00 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $PREFIX_OUT.sbatch.log --wrap " \
python3 /home/jz286/WES_analysis/LDSPEC/CLI_ldspec.py\
    --job regress\
    --pgen_file $PGEN_FILE\
    --annot_file $ANNOT_FILE\
    --score_file $SCORE_FILE\
    --sumstats_file $SUMSTATS_FILE\
    --avgr_file $AVGR_FILE\
    --prefix_out $PREFIX_OUT\
    --flag_nofil_snp True"
done