# #!/bin/bash

PGEN_FILE=/n/scratch3/users/j/jz286/imp_geno_chimp/ukb_imp_chr1_v3_chimp
AVGR_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot/main.avgr

SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.100123/causal_h2g50_p20
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.100123/causal_h2g50_p10
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.100123/causal_h2g20_p20
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.100123/causal_pos_h2g50_p20
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.100123/causal_neg_h2g50_p20
SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.100123/null_h2g50_p20
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.100123/null_h2g50_p10
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.100123/null_h2g20_p20


ANALYSIS=bsl
ANNOT_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001/baseline_annot/baseline_165annots_chr@.annot.gz
SCORE_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/baseline.c@_score.tsv.gz

# ANALYSIS=prox_gene_fct_all_ld
# ANNOT_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.analysis_imp_geno_chimp/reg_annot_file/reg_annot_file.prox_gene_fct_all_ld.txt
# SCORE_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/baseline.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/pannot_prox_ld.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/pannot_gene_ld.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/pannot_baseline_0_100_ld.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/pannot_baseline_0_1000_ld.c@_score.tsv.gz


OUT_PATH=${SIMU_PATH}.${ANALYSIS}

[ -d $OUT_PATH ] || mkdir $OUT_PATH

for i_rep in {1..49}
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