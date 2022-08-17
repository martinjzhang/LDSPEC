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

# Gencode : gene
for CHROM in {2,3,16}
# for CHROM in 16
do
PGEN_FILE=/n/scratch3/users/j/jz286/imp_geno/ukb_imp_chr${CHROM}_v3
OUT_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/pannot/gene
# python get_pannot_s1_gencode.py --pgen_file $PGEN_FILE --pannot gene --out_path $OUT_PATH
sbatch -p short -t 0-03:00 -n 1 -c 1 --mem=128000 --open-mode=truncate -o $OUT_PATH/gene.chr${CHROM}.sbatch.log --wrap " \
python get_pannot_s3_gencode.py --pgen_file $PGEN_FILE --pannot gene --out_path $OUT_PATH"
done

# # Gencode : exon
# for CHROM in {1..21}
# # for CHROM in 22
# do
# PGEN_FILE=/n/scratch3/users/j/jz286/imp_geno/ukb_imp_chr${CHROM}_v3
# OUT_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/pannot_gencode
# # python get_pannot_s1_gene.py --pgen_file $PGEN_FILE --pannot exon --out_path $OUT_PATH
# sbatch -p short -t 0-01:00 -n 1 -c 1 --mem=32000 --open-mode=truncate -o $OUT_PATH/exon.chr${CHROM}.sbatch.log --wrap " \
# python get_pannot_s1_gencode.py --pgen_file $PGEN_FILE --pannot exon --out_path $OUT_PATH"
# done

# # Gencode : exonic_gene
# for CHROM in {1..21}
# # for CHROM in 22
# do
# PGEN_FILE=/n/scratch3/users/j/jz286/imp_geno/ukb_imp_chr${CHROM}_v3
# OUT_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/pannot_gencode
# # python get_pannot_s1_gene.py --pgen_file $PGEN_FILE --pannot exonic_gene --out_path $OUT_PATH
# sbatch -p short -t 0-01:00 -n 1 -c 1 --mem=32000 --open-mode=truncate -o $OUT_PATH/exonic_gene.chr${CHROM}.sbatch.log --wrap " \
# python get_pannot_s1_gencode.py --pgen_file $PGEN_FILE --pannot exonic_gene --out_path $OUT_PATH"
# done

# # Gencode : protein_domain
# for CHROM in {1..22}
# # for CHROM in 22
# do
# PGEN_FILE=/n/scratch3/users/j/jz286/imp_geno/ukb_imp_chr${CHROM}_v3
# OUT_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/pannot_gencode
# # python get_pannot_s1_gene.py --pgen_file $PGEN_FILE --pannot exonic_gene --out_path $OUT_PATH
# sbatch -p short -t 0-01:00 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $OUT_PATH/protein_domain.chr${CHROM}.sbatch.log --wrap " \
# python get_pannot_s2_gencode.py --pgen_file $PGEN_FILE --pannot protein_domain --out_path $OUT_PATH"
# done

# # cS2G_all
# for CHROM in {1..22}
# do
# PGEN_FILE=/n/scratch3/users/j/jz286/imp_geno/ukb_imp_chr${CHROM}_v3
# OUT_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/pannot_cS2G
# # python get_pannot_s3_cS2G.py --pgen_file $PGEN_FILE --pannot cS2G_all --out_path $OUT_PATH
# sbatch -p short -t 0-00:30 -n 1 -c 1 --mem=32000 --open-mode=truncate -o $OUT_PATH/cS2G_all.chr${CHROM}.sbatch.log --wrap "python get_pannot_s3_cS2G.py --pgen_file $PGEN_FILE --pannot cS2G_all --out_path $OUT_PATH"
# done

# # cS2G_promoter
# for CHROM in {1..22}
# do
# PGEN_FILE=/n/scratch3/users/j/jz286/imp_geno/ukb_imp_chr${CHROM}_v3
# OUT_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/pannot_cS2G
# # python get_pannot_s3_cS2G.py --pgen_file $PGEN_FILE --pannot cS2G_promoter --out_path $OUT_PATH
# sbatch -p short -t 0-00:30 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $OUT_PATH/cS2G_promoter.chr${CHROM}.sbatch.log --wrap "python get_pannot_s3_cS2G.py --pgen_file $PGEN_FILE --pannot cS2G_promoter --out_path $OUT_PATH"
# done

# # cS2G_other
# for CHROM in {1..22}
# do
# PGEN_FILE=/n/scratch3/users/j/jz286/imp_geno/ukb_imp_chr${CHROM}_v3
# OUT_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/pannot_cS2G
# # python get_pannot_s3_cS2G.py --pgen_file $PGEN_FILE --pannot cS2G_other --out_path $OUT_PATH
# sbatch -p short -t 0-00:30 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $OUT_PATH/cS2G_other.chr${CHROM}.sbatch.log --wrap "python get_pannot_s3_cS2G.py --pgen_file $PGEN_FILE --pannot cS2G_other --out_path $OUT_PATH"
# done

# # Proxy
# for CHROM in {1..22}
# # for CHROM in 21
# do
# PGEN_FILE=/n/scratch3/users/j/jz286/imp_geno/ukb_imp_chr${CHROM}_v3
# OUT_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/pannot_proxy_1000_10000
# sbatch -p short -t 0-01:00 -n 1 -c 1 --mem=8000 --open-mode=truncate -o $OUT_PATH.proxy.sbatch.log --wrap " \
# python pannot_s2_proxy.py --pgen_file $PGEN_FILE --lb 1000 --ub 10000 --out_path $OUT_PATH"
# done

# # ldp5_proxy_10000
# for CHROM in {1..22}
# do
# PGEN_FILE=/n/scratch3/users/j/jz286/imp_geno/ukb_imp_chr${CHROM}_v3
# LD_FILE=/n/scratch3/users/j/jz286/imp_geno.gdreg_ld/ukb_imp_v3.@_ld.npz
# OUT_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/pannot_ldp5_proxy_10000
# # python pannot_s3_ld.py --pgen_file $PGEN_FILE --ld_file $LD_FILE --out_path $OUT_PATH
# sbatch -p short -t 0-6:00 -n 1 -c 1 --mem=64000 --open-mode=truncate -o $OUT_PATH/ldp5_proxy_10000.chr${CHROM}.sbatch.log --wrap " \
# python pannot_s3_ld.py --pgen_file $PGEN_FILE --ld_file $LD_FILE --out_path $OUT_PATH"
# done