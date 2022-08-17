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

for CHROM in {1..21}
# for CHROM in 22
do
vcf_file=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/vcf/ukb_imp_chr${CHROM}_v3.vcf
vcf_out=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/baseline_annot/vep/ukb_imp_chr${CHROM}_v3.vep_annot
dir_cache=/n/groups/price/martin/vep

# /home/jz286/code/ensembl-vep/vep -i $vcf_file --format vcf\
#     -o $vcf_out --force_overwrite --tab\
#     -a GRCh37 --cache --dir_cache $dir_cache\
#     --sift b --polyphen b --domains --transcript_version --protein --symbol
    
sbatch -p short -t 0-02:00 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $vcf_out.sbatch.log --wrap "\
/home/jz286/code/ensembl-vep/vep -i $vcf_file --format vcf\
    -o $vcf_out --force_overwrite --tab\
    -a GRCh37 --cache --dir_cache $dir_cache\
    --sift b --polyphen b --domains --transcript_version --protein --symbol\
"
done