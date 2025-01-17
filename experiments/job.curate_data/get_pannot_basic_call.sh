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

for CHROM in {1..22}
# for CHROM in 22
do
PGEN_FILE=/n/scratch3/users/j/jz286/imp_geno_chimp/ukb_imp_chr${CHROM}_v3_chimp
LD_FILE=/n/scratch3/users/j/jz286/imp_geno_chimp.gdreg_ld_1e7/ukb_imp_v3_chimp.@_ld.npz
SNP_RANGE_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ukb_imp_v3.snp_range.txt
MAF_BIN_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.curate_data/maf_bin_file.basic.txt

# # prox_0_100, common & lf
# OUT_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot_basic/prox_0_100
# sbatch -p short -t 0-03:00 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $OUT_PATH/chr${CHROM}.sbatch.log --wrap " \
# python3 /home/jz286/WES_analysis/LDSPEC/CLI_pannot_basic.py \
#     --pgen_file $PGEN_FILE\
#     --dist_lb 0 --dist_ub 100\
#     --ld_lb -1 --ld_ub 1\
#     --maf_bin_file $MAF_BIN_FILE\
#     --out_path $OUT_PATH"
    
# # prox_0_100, common & lf, pos-LD
# OUT_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot_basic/prox_0_100_ld_pos
# sbatch -p short -t 0-06:00 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $OUT_PATH/chr${CHROM}.sbatch.log --wrap " \
# python3 /home/jz286/WES_analysis/LDSPEC/CLI_pannot_basic.py \
#     --pgen_file $PGEN_FILE\
#     --ld_file $LD_FILE\
#     --snp_range_file $SNP_RANGE_FILE\
#     --dist_lb 0 --dist_ub 100\
#     --ld_lb 0 --ld_ub 1\
#     --maf_bin_file $MAF_BIN_FILE\
#     --out_path $OUT_PATH"
    
# # prox_0_100, common & lf, neg-LD
# OUT_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot_basic/prox_0_100_ld_neg
# sbatch -p short -t 0-06:00 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $OUT_PATH/chr${CHROM}.sbatch.log --wrap " \
# python3 /home/jz286/WES_analysis/LDSPEC/CLI_pannot_basic.py \
#     --pgen_file $PGEN_FILE\
#     --ld_file $LD_FILE\
#     --snp_range_file $SNP_RANGE_FILE\
#     --dist_lb 0 --dist_ub 100\
#     --ld_lb -1 --ld_ub 0\
#     --maf_bin_file $MAF_BIN_FILE\
#     --out_path $OUT_PATH"
    
# # prox_0_1000, common & lf
# OUT_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot_basic/prox_0_1000
# sbatch -p short -t 0-03:00 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $OUT_PATH/chr${CHROM}.sbatch.log --wrap " \
# python3 /home/jz286/WES_analysis/LDSPEC/CLI_pannot_basic.py \
#     --pgen_file $PGEN_FILE\
#     --dist_lb 0 --dist_ub 1000\
#     --ld_lb -1 --ld_ub 1\
#     --maf_bin_file $MAF_BIN_FILE\
#     --out_path $OUT_PATH"
    
# # prox_0_1000, common & lf, pos-LD
# OUT_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot_basic/prox_0_1000_ld_pos
# sbatch -p short -t 0-06:00 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $OUT_PATH/chr${CHROM}.sbatch.log --wrap " \
# python3 /home/jz286/WES_analysis/LDSPEC/CLI_pannot_basic.py \
#     --pgen_file $PGEN_FILE\
#     --ld_file $LD_FILE\
#     --snp_range_file $SNP_RANGE_FILE\
#     --dist_lb 0 --dist_ub 1000\
#     --ld_lb 0 --ld_ub 1\
#     --maf_bin_file $MAF_BIN_FILE\
#     --out_path $OUT_PATH"
    
# # prox_0_1000, common & lf, neg-LD
# OUT_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot_basic/prox_0_1000_ld_neg
# sbatch -p short -t 0-06:00 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $OUT_PATH/chr${CHROM}.sbatch.log --wrap " \
# python3 /home/jz286/WES_analysis/LDSPEC/CLI_pannot_basic.py \
#     --pgen_file $PGEN_FILE\
#     --ld_file $LD_FILE\
#     --snp_range_file $SNP_RANGE_FILE\
#     --dist_lb 0 --dist_ub 1000\
#     --ld_lb -1 --ld_ub 0\
#     --maf_bin_file $MAF_BIN_FILE\
#     --out_path $OUT_PATH"
    
# prox_100_1000, common & lf
OUT_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot_basic/prox_100_1000
sbatch -p short -t 0-03:00 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $OUT_PATH/chr${CHROM}.sbatch.log --wrap " \
python3 /home/jz286/WES_analysis/LDSPEC/CLI_pannot_basic.py \
    --pgen_file $PGEN_FILE\
    --dist_lb 100 --dist_ub 1000\
    --ld_lb -1 --ld_ub 1\
    --maf_bin_file $MAF_BIN_FILE\
    --out_path $OUT_PATH"
    
# prox_100_1000, common & lf, pos-LD
OUT_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot_basic/prox_100_1000_ld_pos
sbatch -p short -t 0-06:00 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $OUT_PATH/chr${CHROM}.sbatch.log --wrap " \
python3 /home/jz286/WES_analysis/LDSPEC/CLI_pannot_basic.py \
    --pgen_file $PGEN_FILE\
    --ld_file $LD_FILE\
    --snp_range_file $SNP_RANGE_FILE\
    --dist_lb 100 --dist_ub 1000\
    --ld_lb 0 --ld_ub 1\
    --maf_bin_file $MAF_BIN_FILE\
    --out_path $OUT_PATH"
    
# prox_100_1000, common & lf, neg-LD
OUT_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot_basic/prox_100_1000_ld_neg
sbatch -p short -t 0-06:00 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $OUT_PATH/chr${CHROM}.sbatch.log --wrap " \
python3 /home/jz286/WES_analysis/LDSPEC/CLI_pannot_basic.py \
    --pgen_file $PGEN_FILE\
    --ld_file $LD_FILE\
    --snp_range_file $SNP_RANGE_FILE\
    --dist_lb 100 --dist_ub 1000\
    --ld_lb -1 --ld_ub 0\
    --maf_bin_file $MAF_BIN_FILE\
    --out_path $OUT_PATH"
    
# # prox_0_10000, common & lf
# OUT_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot_basic/prox_0_10000
# sbatch -p short -t 0-03:00 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $OUT_PATH/chr${CHROM}.sbatch.log --wrap " \
# python3 /home/jz286/WES_analysis/LDSPEC/CLI_pannot_basic.py \
#     --pgen_file $PGEN_FILE\
#     --dist_lb 0 --dist_ub 10000\
#     --ld_lb -1 --ld_ub 1\
#     --maf_bin_file $MAF_BIN_FILE\
#     --out_path $OUT_PATH"
    
# # prox_0_10000, common & lf, pos-LD
# OUT_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot_basic/prox_0_10000_ld_pos
# sbatch -p short -t 0-06:00 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $OUT_PATH/chr${CHROM}.sbatch.log --wrap " \
# python3 /home/jz286/WES_analysis/LDSPEC/CLI_pannot_basic.py \
#     --pgen_file $PGEN_FILE\
#     --ld_file $LD_FILE\
#     --snp_range_file $SNP_RANGE_FILE\
#     --dist_lb 0 --dist_ub 10000\
#     --ld_lb 0 --ld_ub 1\
#     --maf_bin_file $MAF_BIN_FILE\
#     --out_path $OUT_PATH"
    
# # prox_0_10000, common & lf, neg-LD
# OUT_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot_basic/prox_0_10000_ld_neg
# sbatch -p short -t 0-06:00 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $OUT_PATH/chr${CHROM}.sbatch.log --wrap " \
# python3 /home/jz286/WES_analysis/LDSPEC/CLI_pannot_basic.py \
#     --pgen_file $PGEN_FILE\
#     --ld_file $LD_FILE\
#     --snp_range_file $SNP_RANGE_FILE\
#     --dist_lb 0 --dist_ub 10000\
#     --ld_lb -1 --ld_ub 0\
#     --maf_bin_file $MAF_BIN_FILE\
#     --out_path $OUT_PATH"
    
# prox_1000_10000, common & lf
OUT_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot_basic/prox_1000_10000
sbatch -p short -t 0-03:00 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $OUT_PATH/chr${CHROM}.sbatch.log --wrap " \
python3 /home/jz286/WES_analysis/LDSPEC/CLI_pannot_basic.py \
    --pgen_file $PGEN_FILE\
    --dist_lb 1000 --dist_ub 10000\
    --ld_lb -1 --ld_ub 1\
    --maf_bin_file $MAF_BIN_FILE\
    --out_path $OUT_PATH"
    
# prox_1000_10000, common & lf, pos-LD
OUT_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot_basic/prox_1000_10000_ld_pos
sbatch -p short -t 0-06:00 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $OUT_PATH/chr${CHROM}.sbatch.log --wrap " \
python3 /home/jz286/WES_analysis/LDSPEC/CLI_pannot_basic.py \
    --pgen_file $PGEN_FILE\
    --ld_file $LD_FILE\
    --snp_range_file $SNP_RANGE_FILE\
    --dist_lb 1000 --dist_ub 10000\
    --ld_lb 0 --ld_ub 1\
    --maf_bin_file $MAF_BIN_FILE\
    --out_path $OUT_PATH"
    
# prox_1000_10000, common & lf, neg-LD
OUT_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot_basic/prox_1000_10000_ld_neg
sbatch -p short -t 0-06:00 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $OUT_PATH/chr${CHROM}.sbatch.log --wrap " \
python3 /home/jz286/WES_analysis/LDSPEC/CLI_pannot_basic.py \
    --pgen_file $PGEN_FILE\
    --ld_file $LD_FILE\
    --snp_range_file $SNP_RANGE_FILE\
    --dist_lb 1000 --dist_ub 10000\
    --ld_lb -1 --ld_ub 0\
    --maf_bin_file $MAF_BIN_FILE\
    --out_path $OUT_PATH"
done