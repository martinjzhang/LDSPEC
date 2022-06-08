#!/bin/bash

PGEN_FILE=/n/groups/price/martin/WES_analysis/toy_1K/chr@_v1.SPB.hg19.toy_1K
ANNOT_FILE=/n/groups/price/martin/WES_analysis/toy_1K/toy.annot.gz,/n/groups/price/martin/WES_analysis/toy_1K/toy.pannot.gz,/n/groups/price/martin/WES_analysis/toy_1K/toy.pannot_hr.gz



for REP in {0..19}
do
SUMSTATS_FILE=/n/groups/price/martin/WES_analysis/toy_1K/sanity_rep${REP}.sumstats.gz
PREFIX_OUT=/n/groups/price/martin/WES_analysis/toy_1K/results/reg/sanity_rep${REP}

sbatch -p short -t 0-00:20 -n 1 -c 1 --mem=8000 --open-mode=truncate -e $PREFIX_OUT.sbatch.err -o $PREFIX_OUT.sbatch.log --wrap " \
python3 /home/jz286/WES_analysis/GDReg/run_gdreg.py\
    --job regress\
    --pgen_file $PGEN_FILE\
    --ld_file /n/groups/price/martin/WES_analysis/toy_1K/results/top_1K_chr@.ld.npz\
    --sumstats_file $SUMSTATS_FILE\
    --annot_file $ANNOT_FILE\
    --random_seed 0\
    --memory 2048\
    --prefix_out $PREFIX_OUT"
done