#!/bin/bash
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-06:00
#SBATCH --array=2-22
#SBATCH -p short
#SBATCH --mem=32000
#SBATCH -o /home/jz286/WES_analysis/GDReg/experiments/job_info/hostname_%j.out  
#SBATCH -e /home/jz286/WES_analysis/GDReg/experiments/job_info/hostname_%j.err 
#SBATCH --mail-type=NONE#SBATCH --mail-type=NONE

# # Step1: Create the position file for tabix: I just copied the python files here
# import pandas as pd
# import numpy as np
    
# DATA_PATH='/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/baseline_annot/bim'
# OUT_PATH='/n/scratch3/users/j/jz286/temp'
# for CHROM in np.arange(1,23):
#     df_bim=pd.read_csv(DATA_PATH+'/ukb_imp_chr%d_v3.bim'%CHROM,
#                        delim_whitespace=True,header=None,index_col=False)
#     df_bim.columns=['CHR','SNP','CM','BP','A1','A2']
#     df_bim = df_bim[['CHR', 'BP']]
#     df_bim.to_csv(OUT_PATH+'/ukb_imp_chr%d_v3.bim'%CHROM, sep='\t', header=False, index=False)



# # Step2: Get the CADD annotion for the region using tabix
# # Each file takes ~3hrs to complete
# DATA_PATH=/n/scratch3/users/j/jz286/CADD
# CHROM=$SLURM_ARRAY_TASK_ID
# POS_FILE=/n/scratch3/users/j/jz286/temp/ukb_imp_chr${CHROM}_v3.bim
# OUT_FILE=/n/scratch3/users/j/jz286/temp/ukb_imp_chr${CHROM}_v3.annot
# tabix -h -R $POS_FILE $DATA_PATH/whole_genome_SNVs_inclAnno.tsv.gz > $OUT_FILE



# Step3: Use get_CADD_annot.py to get the final (allele-specific) .annot file

for CHROM in {1..22}
# for CHROM in 21
do
BIM_FILE=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/baseline_annot/bim/ukb_imp_chr${CHROM}_v3.bim
CADD_FILE=/n/scratch3/users/j/jz286/temp/ukb_imp_chr${CHROM}_v3.annot
OUT_FILE=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/baseline_annot/CADD/ukb_imp_chr${CHROM}_v3.CADD
python baseline_s4_CADD_annot.py --bim_file $BIM_FILE --cadd_file $CADD_FILE --out_file $OUT_FILE
done

