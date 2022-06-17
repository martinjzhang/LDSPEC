#!/bin/bash

PGEN_FILE=/n/groups/price/martin/WES_analysis/toy_10K/chr@_v1.SPB.hg19.toy_10K
ANNOT_FILE=/n/groups/price/martin/WES_analysis/toy_10K/toy.annot.gz,/n/groups/price/martin/WES_analysis/toy_10K/toy.pannot.gz,/n/groups/price/martin/WES_analysis/toy_10K/toy.pannot_hr.gz



for REP in {0..19}
# for REP in 0
do
SCORE_FILE=/n/groups/price/martin/WES_analysis/toy_10K/gdreg_file/toy_10K.@_score.tsv.gz
SUMSTATS_FILE=/n/groups/price/martin/WES_analysis/toy_10K/sanity_rep${REP}.sumstats.gz
PREFIX_OUT=/n/groups/price/martin/WES_analysis/toy_10K/gdreg_res/sanity_rep${REP}

# python3 /home/jz286/WES_analysis/GDReg/run_gdreg.py\
#     --job regress\
#     --pgen_file $PGEN_FILE\
#     --score_file $SCORE_FILE\
#     --sumstats_file $SUMSTATS_FILE\
#     --annot_file $ANNOT_FILE\
#     --random_seed 0\
#     --memory 2048\
#     --prefix_out $PREFIX_OUT

sbatch -p short -t 0-00:10 -n 1 -c 1 --mem=4000 --open-mode=truncate -e $PREFIX_OUT.sbatch.err -o $PREFIX_OUT.sbatch.log --wrap " \
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