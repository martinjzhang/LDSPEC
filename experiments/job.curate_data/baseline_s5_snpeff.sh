#!/bin/bash
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-02:00
#SBATCH --array=1-22
#SBATCH -p short
#SBATCH --mem=16000
#SBATCH -o /home/jz286/WES_analysis/wes_rare/experiments/job_info/hostname_%j.out  
#SBATCH -e /home/jz286/WES_analysis/wes_rare/experiments/job_info/hostname_%j.err 
#SBATCH --mail-type=NONE#SBATCH --mail-type=NONE

CHROM=$SLURM_ARRAY_TASK_ID
# CHROM=22


for CHROM in {1..20}
# for CHROM in 22
do
VCF_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/vcf
vcf_file=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/vcf/ukb_imp_chr${CHROM}_v3.vcf
vcf_out=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/baseline_annot/snfeff/ukb_imp_chr${CHROM}_v3.vcf
sbatch -p short -t 0-01:00 -n 1 -c 1 --mem=8000 --open-mode=truncate -o $PREFIX_OUT.convert_bed.sbatch.log --wrap " \
java -Xmx4g -jar /home/jz286/snpEff/snpEff.jar -c /home/jz286/snpEff/snpEff.config -lof -v GRCh37.75 $vcf_file > $vcf_out"
done