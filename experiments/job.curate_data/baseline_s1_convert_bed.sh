#!/bin/bash
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-01:00
#SBATCH --array=1
#SBATCH -p short
#SBATCH --mem=32000
#SBATCH -o /home/jz286/WES_analysis/wes_rare/experiments/job_info/hostname_%j.out  
#SBATCH -e /home/jz286/WES_analysis/wes_rare/experiments/job_info/hostname_%j.err 
#SBATCH --mail-type=NONE#SBATCH --mail-type=NONE
# CHROM=22

for CHROM in {1..22}
# for CHROM in 22
do

# echo $CHROM
PGEN_FILE=/n/scratch3/users/j/jz286/imp_geno/ukb_imp_chr${CHROM}_v3
PREFIX_OUT=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/baseline_annot/bed/ukb_imp_chr${CHROM}_v3

sbatch -p short -t 0-1:00 -n 1 -c 1 --mem=32000 --open-mode=truncate -o $PREFIX_OUT.convert_bed.sbatch.log --wrap " \
python baseline_s1_convert_bed.py --pgen_file $PGEN_FILE --prefix_out $PREFIX_OUT"
    
# python baseline_s1_convert_bed.py --pgen_file $PGEN_FILE --prefix_out $PREFIX_OUT

done