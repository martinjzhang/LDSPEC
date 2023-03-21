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

# # Gencode : gene
# for CHROM in {2,3,16}
# # for CHROM in 16
# do
# PGEN_FILE=/n/scratch3/users/j/jz286/imp_geno/ukb_imp_chr${CHROM}_v3
# OUT_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/pannot/gene
# # python get_pannot_s1_gencode.py --pgen_file $PGEN_FILE --pannot gene --out_path $OUT_PATH
# sbatch -p short -t 0-03:00 -n 1 -c 1 --mem=128000 --open-mode=truncate -o $OUT_PATH/gene.chr${CHROM}.sbatch.log --wrap " \
# python get_pannot_s3_gencode.py --pgen_file $PGEN_FILE --pannot gene --out_path $OUT_PATH"
# done

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
# # for CHROM in 22
# do
# PGEN_FILE=/n/scratch3/users/j/jz286/imp_geno/ukb_imp_chr${CHROM}_v3
# # OUT_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/pannot/proxy_0_1000
# OUT_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/pannot/proxy_0_10000
# # python get_pannot_s1_proxy.py --pgen_file $PGEN_FILE --lb 0 --ub 1000 --out_path $OUT_PATH
# sbatch -p short -t 0-01:00 -n 1 -c 1 --mem=8000 --open-mode=truncate -o $OUT_PATH.proxy.sbatch.log --wrap " \
# python get_pannot_s1_proxy.py --pgen_file $PGEN_FILE --lb 0 --ub 10000 --out_path $OUT_PATH"
# done

# ldp5_proxy
for CHROM in {1..21}
# for CHROM in 22
do
PGEN_FILE=/n/scratch3/users/j/jz286/imp_geno/ukb_imp_chr${CHROM}_v3
LD_FILE=/n/scratch3/users/j/jz286/imp_geno.gdreg_ld/ukb_imp_v3.@_ld.npz
OUT_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/pannot/ldp5_proxy_0_1000
# OUT_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/pannot/pannot_ldp5_proxy_100_1000
# OUT_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/pannot/pannot_ldp5_proxy_1000_10000
# python get_pannot_s2_ld.py --pgen_file $PGEN_FILE --lb 0 --ub 100 --ld_file $LD_FILE --out_path $OUT_PATH
sbatch -p short -t 0-3:00 -n 1 -c 1 --mem=32000 --open-mode=truncate -o $OUT_PATH/chr${CHROM}.sbatch.log --wrap " \
python get_pannot_s2_ld.py --pgen_file $PGEN_FILE --lb 0 --ub 1000 --ld_file $LD_FILE --out_path $OUT_PATH"
done

# # Loop
# # for CHROM in {1..21}
# for CHROM in 22
# do
# PGEN_FILE=/n/scratch3/users/j/jz286/imp_geno/ukb_imp_chr${CHROM}_v3
# LOOP_FILE=/n/groups/price/martin/data_GDREG/gene_annotation/3dgenome_loops/Dixon_2015.H1-ESC.hg19.peakachu-merged.loops
# OUT_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/pannot/Dixon_2015_H1_ESC
# LOOP_FILE=/n/groups/price/martin/data_GDREG/gene_annotation/3dgenome_loops/Dixon_2015.H1-MES.hg19.peakachu-merged.loops
# OUT_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/pannot/Dixon_2015_H1_MES
# # LOOP_FILE=/n/groups/price/martin/data_GDREG/gene_annotation/3dgenome_loops/Dixon_2015.H1-MSC.hg19.peakachu-merged.loops
# # OUT_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/pannot/Dixon_2015_H1_MSC
# # LOOP_FILE=/n/groups/price/martin/data_GDREG/gene_annotation/3dgenome_loops/Dixon_2015.H1-NPC.hg19.peakachu-merged.loops
# # OUT_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/pannot/Dixon_2015_H1_NPC
# # LOOP_FILE=/n/groups/price/martin/data_GDREG/gene_annotation/3dgenome_loops/Dixon_2015.H1-TRO.hg19.peakachu-merged.loops
# # OUT_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/pannot/Dixon_2015_H1_TRO
# # python get_pannot_s5_loop.py --pgen_file $PGEN_FILE --loop_file $LOOP_FILE --out_path $OUT_PATH
# sbatch -p short -t 0-1:00 -n 1 -c 1 --mem=32000 --open-mode=truncate -o $OUT_PATH/chr${CHROM}.sbatch.log --wrap " \
# python get_pannot_s5_loop.py --pgen_file $PGEN_FILE --loop_file $LOOP_FILE --out_path $OUT_PATH"
# done

# # Baseline
# for CHROM in {1..22}
# # for CHROM in 22
# do
# PGEN_FILE=/n/scratch3/users/j/jz286/imp_geno/ukb_imp_chr${CHROM}_v3
# ANNOT_FILE=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/baseline_annot/baseline_165annots_chr${CHROM}.annot.gz
# OUT_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/pannot/baseline_0_100
# OUT_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/pannot/baseline_100_1000
# OUT_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/pannot/baseline_1000_10000
# OUT_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/pannot/baseline_0_1000
# OUT_PATH=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/pannot/baseline_0_10000
# # python get_pannot_s6_baseline.py --pgen_file $PGEN_FILE --annot_file $ANNOT_FILE --lb 0 --ub 100 --out_path $OUT_PATH
# sbatch -p short -t 0-03:00 -n 1 -c 1 --mem=32000 --open-mode=truncate -o $OUT_PATH/chr${CHROM}.sbatch.log --wrap " \
# python get_pannot_s6_baseline.py --pgen_file $PGEN_FILE --annot_file $ANNOT_FILE --lb 0 --ub 10000 --out_path $OUT_PATH"
# done