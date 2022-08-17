#!/bin/bash
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-02:00
#SBATCH --array= 1
#SBATCH -p short
#SBATCH --mem=16000
#SBATCH -o /home/jz286/WES_analysis/wes_rare/experiments/job_info/hostname_%j.out  
#SBATCH -e /home/jz286/WES_analysis/wes_rare/experiments/job_info/hostname_%j.err 
#SBATCH --mail-type=NONE#SBATCH --mail-type=NONE

# CHROM=$SLURM_ARRAY_TASK_ID
# CHROM=22

for CHROM in {10..11}
# for CHROM in 22
do
# echo $CHROM
INPUT_FILE=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/baseline_annot/bim/ukb_imp_chr${CHROM}_v3.bim
OUTPUT_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/baseline_annot/nucleotide_diversity
sbatch -p short -t 0-3:00 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $PREFIX_OUT.nucleotide_diversity.sbatch.log --wrap " \
Rscript baseline_s2_nucleotide_diversity.R $CHROM $INPUT_FILE $OUTPUT_PATH"
done