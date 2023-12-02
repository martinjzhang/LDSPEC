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


for DIST_UB in {100,1000}
do
for CHROM in {1..22}
# for CHROM in 22
do
PGEN_FILE=/n/scratch3/users/j/jz286/imp_geno_chimp/ukb_imp_chr${CHROM}_v3_chimp
ANNOT_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001/baseline_annot/baseline_165annots_chr${CHROM}.annot.gz
LD_FILE=/n/scratch3/users/j/jz286/imp_geno_chimp.gdreg_ld_1e7/ukb_imp_v3_chimp.@_ld.npz
SNP_RANGE_FILE=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ukb_imp_v3.snp_range.txt
MAF_BIN_FILE=/home/jz286/WES_analysis/LDSPEC/experiments/job.curate_data/maf_bin_file.basic.txt

# common & lf
OUT_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot_fct/baseline_0_${DIST_UB}
[ -d $OUT_PATH ] || mkdir $OUT_PATH
sbatch -p short -t 0-03:00 -n 1 -c 1 --mem=32000 --open-mode=truncate -o $OUT_PATH/chr${CHROM}.sbatch.log --wrap " \
python3 /home/jz286/WES_analysis/LDSPEC/CLI_pannot_fct.py \
    --pgen_file $PGEN_FILE\
    --annot_file $ANNOT_FILE\
    --dist_lb 0 --dist_ub ${DIST_UB}\
    --maf_bin_file $MAF_BIN_FILE\
    --out_path $OUT_PATH"
    
# common & lf, pos-LD
OUT_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot_fct/baseline_0_${DIST_UB}_ld_pos
[ -d $OUT_PATH ] || mkdir $OUT_PATH
sbatch -p short -t 0-12:00 -n 1 -c 1 --mem=64000 --open-mode=truncate -o $OUT_PATH/chr${CHROM}.sbatch.log --wrap " \
python3 /home/jz286/WES_analysis/LDSPEC/CLI_pannot_fct.py \
    --pgen_file $PGEN_FILE\
    --annot_file $ANNOT_FILE\
    --ld_file $LD_FILE\
    --snp_range_file $SNP_RANGE_FILE\
    --dist_lb 0 --dist_ub ${DIST_UB}\
    --ld_lb 0 --ld_ub 2\
    --maf_bin_file $MAF_BIN_FILE\
    --out_path $OUT_PATH"
    
# common & lf, neg-LD
OUT_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/pannot_fct/baseline_0_${DIST_UB}_ld_neg
[ -d $OUT_PATH ] || mkdir $OUT_PATH
sbatch -p short -t 0-12:00 -n 1 -c 1 --mem=64000 --open-mode=truncate -o $OUT_PATH/chr${CHROM}.sbatch.log --wrap " \
python3 /home/jz286/WES_analysis/LDSPEC/CLI_pannot_fct.py \
    --pgen_file $PGEN_FILE\
    --annot_file $ANNOT_FILE\
    --ld_file $LD_FILE\
    --snp_range_file $SNP_RANGE_FILE\
    --dist_lb 0 --dist_ub ${DIST_UB}\
    --ld_lb -2 --ld_ub 0\
    --maf_bin_file $MAF_BIN_FILE\
    --out_path $OUT_PATH"
done
done