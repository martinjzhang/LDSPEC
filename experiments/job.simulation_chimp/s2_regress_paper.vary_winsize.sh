# #!/bin/bash

PGEN_FILE=/n/scratch3/users/j/jz286/imp_geno_chimp/ukb_imp_chr1_v3_chimp
AVGR_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot/main.avgr

SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.100123/null_h2g50_p20
# SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.100123/causal_h2g50_p20

ANALYSIS=prox_gene_fct_all_ld_1e6
ANNOT_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.analysis_imp_geno_chimp/reg_annot_file/reg_annot_file.prox_gene_fct_all_ld.txt
SCORE_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/chr1_vary_winsize/baseline.1e6.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/chr1_vary_winsize/pannot_prox_ld.1e6.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/chr1_vary_winsize/pannot_gene_ld.1e6.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/chr1_vary_winsize/pannot_baseline_0_100_ld.1e6.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/chr1_vary_winsize/pannot_baseline_0_1000_ld.1e6.c@_score.tsv.gz

ANALYSIS=prox_gene_fct_all_ld_3e6
ANNOT_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.analysis_imp_geno_chimp/reg_annot_file/reg_annot_file.prox_gene_fct_all_ld.txt
SCORE_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/chr1_vary_winsize/baseline.3e6.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/chr1_vary_winsize/pannot_prox_ld.3e6.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/chr1_vary_winsize/pannot_gene_ld.3e6.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/chr1_vary_winsize/pannot_baseline_0_100_ld.3e6.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/chr1_vary_winsize/pannot_baseline_0_1000_ld.3e6.c@_score.tsv.gz

ANALYSIS=prox_gene_fct_all_ld_5e6
ANNOT_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.analysis_imp_geno_chimp/reg_annot_file/reg_annot_file.prox_gene_fct_all_ld.txt
SCORE_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/chr1_vary_winsize/baseline.5e6.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/chr1_vary_winsize/pannot_prox_ld.5e6.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/chr1_vary_winsize/pannot_gene_ld.5e6.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/chr1_vary_winsize/pannot_baseline_0_100_ld.5e6.c@_score.tsv.gz,/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldspec_score_file/chr1_vary_winsize/pannot_baseline_0_1000_ld.5e6.c@_score.tsv.gz


OUT_PATH=${SIMU_PATH}.${ANALYSIS}

[ -d $OUT_PATH ] || mkdir $OUT_PATH

for i_rep in {0..49}
# for i_rep in 0
do
SUMSTATS_FILE=${SIMU_PATH}/rep${i_rep}.sumstats.gz
PREFIX_OUT=${OUT_PATH}/rep${i_rep}

sbatch -p short -t 0-03:00 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $PREFIX_OUT.sbatch.log --wrap " \
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