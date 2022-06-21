#!/bin/bash

PGEN_FILE=/n/groups/price/martin/WES_analysis/toy_10K/chr@_v1.SPB.hg19.toy_10K
ANNOT_FILE=/n/groups/price/martin/WES_analysis/toy_10K/toy.annot.gz,/n/groups/price/martin/WES_analysis/toy_10K/toy.pannot.gz,/n/groups/price/martin/WES_analysis/toy_10K/toy.pannot_hr.gz



for i_line in {1..180}
# for i_line in 1
do
TRAIT=$( head -n $i_line "/n/groups/price/martin/WES_analysis/toy_10K/trait_simu_debug/trait_list.txt" | tail -1 )
SCORE_FILE=/n/groups/price/martin/WES_analysis/toy_10K/gdreg_file_score_cross_term/toy_10K.@_score.tsv.gz
SUMSTATS_FILE=/n/groups/price/martin/WES_analysis/toy_10K/trait_simu_debug/${TRAIT}.sumstats.gz
PREFIX_OUT=/n/groups/price/martin/WES_analysis/toy_10K/gdreg_res_simu_debug.with_cross/${TRAIT}

# python3 /home/jz286/WES_analysis/GDReg/run_gdreg.py\
#     --job regress\
#     --pgen_file $PGEN_FILE\
#     --score_file $SCORE_FILE\
#     --sumstats_file $SUMSTATS_FILE\
#     --annot_file $ANNOT_FILE\
#     --random_seed 0\
#     --memory 2048\
#     --flag_cross_term True\
#     --prefix_out $PREFIX_OUT

sbatch -p short -t 0-00:05 -n 1 -c 1 --mem=4000 --open-mode=truncate -o $PREFIX_OUT.sbatch.log --wrap " \
python3 /home/jz286/WES_analysis/GDReg/run_gdreg.py\
    --job regress\
    --pgen_file $PGEN_FILE\
    --score_file $SCORE_FILE\
    --sumstats_file $SUMSTATS_FILE\
    --annot_file $ANNOT_FILE\
    --random_seed 0\
    --memory 2048\
    --flag_cross_term True\
    --prefix_out $PREFIX_OUT"
done


# Without cross term
for i_line in {1..180}
# for i_line in 1
do
TRAIT=$( head -n $i_line "/n/groups/price/martin/WES_analysis/toy_10K/trait_simu_debug/trait_list.txt" | tail -1 )
SCORE_FILE=/n/groups/price/martin/WES_analysis/toy_10K/gdreg_file_score_cross_term/toy_10K.@_score.tsv.gz
SUMSTATS_FILE=/n/groups/price/martin/WES_analysis/toy_10K/trait_simu_debug/${TRAIT}.sumstats.gz
PREFIX_OUT=/n/groups/price/martin/WES_analysis/toy_10K/gdreg_res_simu_debug.without_cross/${TRAIT}

# python3 /home/jz286/WES_analysis/GDReg/run_gdreg.py\
#     --job regress\
#     --pgen_file $PGEN_FILE\
#     --score_file $SCORE_FILE\
#     --sumstats_file $SUMSTATS_FILE\
#     --annot_file $ANNOT_FILE\
#     --random_seed 0\
#     --memory 2048\
#     --prefix_out $PREFIX_OUT

sbatch -p short -t 0-00:05 -n 1 -c 1 --mem=4000 --open-mode=truncate -o $PREFIX_OUT.sbatch.log --wrap " \
python3 /home/jz286/WES_analysis/GDReg/run_gdreg.py\
    --job regress\
    --pgen_file $PGEN_FILE\
    --score_file $SCORE_FILE\
    --sumstats_file $SUMSTATS_FILE\
    --annot_file $ANNOT_FILE\
    --random_seed 0\
    --memory 2048\
    --prefix_out $PREFIX_OUT"
done