#!/bin/bash
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-12:00
#SBATCH --array=1-2
#SBATCH -p short
#SBATCH --mem=8000
#SBATCH -o /home/jz286/job_info/hostname_%j.out  
#SBATCH -e /home/jz286/job_info/hostname_%j.err 
#SBATCH --mail-type=NONE#SBATCH --mail-type=NONE

# progress: 1-13500
FILE_PATH=/n/groups/price/martin/LDSPEC_data/res_cojo_rv1
i_line=$SLURM_ARRAY_TASK_ID # 1-18958
# i_line=$(($SLURM_ARRAY_TASK_ID + 10000))
echo $i_line

# for TRAIT in {biochemistry_AlkalinePhosphatase,biochemistry_AspartateAminotransferase,biochemistry_Cholesterol,biochemistry_Creatinine,biochemistry_IGF1,biochemistry_Phosphate,biochemistry_TotalBilirubin,biochemistry_TotalProtein,biochemistry_VitaminD,blood_PLATELET_COUNT,blood_RBC_DISTRIB_WIDTH,blood_RED_COUNT,blood_WHITE_COUNT,bmd_HEEL_TSCOREz,body_BALDING1}
# for TRAIT in {body_BMIz,body_HEIGHTz,body_WHRadjBMIz,bp_DIASTOLICadjMEDz,cov_EDU_YEARS,disease_ALLERGY_ECZEMA_DIAGNOSED,disease_HYPOTHYROIDISM_SELF_REP,lung_FEV1FVCzSMOKE,lung_FVCzSMOKE,mental_NEUROTICISM,other_MORNINGPERSON,repro_MENARCHE_AGE,repro_MENOPAUSE_AGE,repro_NumberChildrenEverBorn_Pooled}
for TRAIT in body_BMIz
do

SUMSTATS_FILE=/n/groups/price/martin/LDSPEC_data/res_cojo_rv1/sumstats/${TRAIT}.sumstats
# 100-kb window
OUT_PATH=/n/groups/price/martin/LDSPEC_data/res_cojo_rv1/res_cojo_100kb/${TRAIT}
GENE_FILE_PATH=/n/groups/price/martin/LDSPEC_data/res_cojo_rv1/gene_annotation/gene_level_snp_100kb

GENE_LIST_FILE=/n/groups/price/martin/LDSPEC_data/res_cojo_rv1/gene_annotation/gene_level_snp_100kb.txt
GENE_LIST_FILE=/n/groups/price/martin/LDSPEC_data/res_cojo_rv1/gene_annotation/gene_level_snp_100kb_uf.txt

GENE_FILE=$( head -n $i_line ${GENE_LIST_FILE} | tail -1 )
arrIN=($GENE_FILE)
CHR=${arrIN[1]}
GENE_FILE=${arrIN[2]}
BFILE=/n/scratch/users/j/jz286/imp_geno_chimp_bed_MAF01/ukb_imp_chr${CHR}_v3_chimp

echo ${GENE_FILE}

[ -d $OUT_PATH ] || mkdir $OUT_PATH


/home/jz286/code/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1\
    --bfile ${BFILE}\
    --chr ${CHR}\
    --extract ${GENE_FILE_PATH}/${GENE_FILE}\
    --maf 0.01\
    --cojo-file ${SUMSTATS_FILE}\
    --cojo-slct\
    --out ${OUT_PATH}/${GENE_FILE}
done