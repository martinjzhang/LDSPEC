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

PGEN_FILE=/n/groups/price/martin/data_GDREG/toy_10K/chr@_v1.SPB.hg19.toy_10K
ANNOT_FILE=/n/groups/price/martin/data_GDREG/toy_10K/toy.chr@.annot.gz,/n/groups/price/martin/data_GDREG/toy_10K/toy.gene.chr@.pannot_mat.npz,/n/groups/price/martin/data_GDREG/toy_10K/toy.proxy.chr@.pannot_mat.npz

LD_FILE=/n/groups/price/martin/data_GDREG/toy_10K/gdreg_file_ld/toy_10K.@_ld.npz    
PREFIX_OUT=/n/groups/price/martin/data_GDREG/toy_10K/pannot

python3 /home/jz286/WES_analysis/GDReg/run_gdreg.py\
    --job compute_avgr\
    --pgen_file $PGEN_FILE\
    --ld_file $LD_FILE\
    --annot_file $ANNOT_FILE\
    --prefix_out $PREFIX_OUT