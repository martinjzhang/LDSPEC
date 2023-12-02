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

# for pannot in {gene,exon,exonic_gene,protein_domain,cS2G_promoter}
for pannot in exon
do
# for CHROM in {1..3}
for CHROM in 2
do
PGEN_FILE=/n/scratch3/users/j/jz286/imp_geno_chimp/ukb_imp_chr${CHROM}_v3_chimp
LD_FILE=/n/scratch3/users/j/jz286/imp_geno_chimp.gdreg_ld_1e7/ukb_imp_v3_chimp.@_ld.npz
SNP_RANGE_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ukb_imp_v3.snp_range.txt
MAF_BIN_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.curate_data/maf_bin_file.basic.txt

# # common & lf
# OUT_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot_gene/${pannot}
# [ -d $OUT_PATH ] || mkdir $OUT_PATH
# sbatch -p short -t 0-12:00 -n 1 -c 1 --mem=96000 --open-mode=truncate -o $OUT_PATH/chr${CHROM}.sbatch.log --wrap " \
# python3 /home/jz286/WES_analysis/LDSPEC/CLI_pannot_gene.py \
#     --pgen_file $PGEN_FILE\
#     --pannot ${pannot}\
#     --ld_lb -1 --ld_ub 1\
#     --maf_bin_file $MAF_BIN_FILE\
#     --out_path $OUT_PATH"
    
# # common & lf, pos-LD
# OUT_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot_gene/${pannot}_ld_pos
# [ -d $OUT_PATH ] || mkdir $OUT_PATH
# sbatch -p medium -t 0-24:00 -n 1 -c 1 --mem=64000 --open-mode=truncate -o $OUT_PATH/chr${CHROM}.sbatch.log --wrap " \
# python3 /home/jz286/WES_analysis/LDSPEC/CLI_pannot_gene.py \
#     --pgen_file $PGEN_FILE\
#     --pannot ${pannot}\
#     --ld_file $LD_FILE\
#     --snp_range_file $SNP_RANGE_FILE\
#     --ld_lb 0 --ld_ub 2\
#     --maf_bin_file $MAF_BIN_FILE\
#     --out_path $OUT_PATH"
    
# common & lf, neg-LD
OUT_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot_gene/${pannot}_ld_neg
[ -d $OUT_PATH ] || mkdir $OUT_PATH
sbatch -p medium -t 2-00:00 -n 1 -c 1 --mem=64000 --open-mode=truncate -o $OUT_PATH/chr${CHROM}.sbatch.log --wrap " \
python3 /home/jz286/WES_analysis/LDSPEC/CLI_pannot_gene.py \
    --pgen_file $PGEN_FILE\
    --pannot ${pannot}\
    --ld_file $LD_FILE\
    --snp_range_file $SNP_RANGE_FILE\
    --ld_lb -2 --ld_ub 0\
    --maf_bin_file $MAF_BIN_FILE\
    --out_path $OUT_PATH"
done
done