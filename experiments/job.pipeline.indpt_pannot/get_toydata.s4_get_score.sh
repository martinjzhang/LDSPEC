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
    
    
# compute_ld
for i_line in {1..9}
# for i_line in 10
do 
PREFIX_OUT=/n/groups/price/martin/data_GDREG/toy_10K/gdreg_file_score.lazy/toy_10K
SNP_RANGE=$( head -n $i_line "/n/groups/price/martin/data_GDREG/toy_10K/toy_10K.snp_range.txt" | tail -1 )
LD_FILE=/n/groups/price/martin/data_GDREG/toy_10K/gdreg_file_ld/toy_10K.${SNP_RANGE}_ld.npz

# python3 /home/jz286/WES_analysis/GDReg/run_gdreg.py\
#     --job compute_score\
#     --pgen_file $PGEN_FILE\
#     --ld_file $LD_FILE\
#     --annot_file $ANNOT_FILE\
#     --memory 2048\
#     --flag_cross_term True\
#     --prefix_out $PREFIX_OUT
    
sbatch -p short -t 0-00:10 -n 1 -c 1 --mem=4000 --open-mode=truncate -o ${PREFIX_OUT}.${SNP_RANGE}.compute_score.sbatch.log --wrap " \
python3 /home/jz286/WES_analysis/GDReg/run_gdreg.py\
    --job compute_score\
    --pgen_file $PGEN_FILE\
    --ld_file $LD_FILE\
    --annot_file $ANNOT_FILE\
    --memory 2048\
    --flag_cross_term True\
    --prefix_out $PREFIX_OUT"
done