#!/bin/bash

# $ldscore = "python /n/groups/price/steven/soft/ldsc/ldsc.py";
# $options = "--overlap-annot --print-coefficients"; # --print-cov --print-delete-vals";

# $freq    = "--frqfile-chr /n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC.";
# $weights = "--w-ld-chr /n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/weights/weights.hm3_noMHC.";
# $annots1 = "/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baselineLD_v2.2/baselineLD.";



TRAIT=rep0
sumstats=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001_chimp/simulation.040123/null_h2g50_p20/${TRAIT}.sumstats.gz
PATH_OUT=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001_chimp/simulation.040123/ldsc
PREFIX_OUT=${PATH_OUT}/${TRAIT}.bslLF165

python /n/groups/price/steven/soft/ldsc/ldsc.py \
    --h2 $sumstats \
    --ref-ld /n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/ldsc_res/ld_file/bslLF165.1 \
    --w-ld /n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/ldsc_res/ld_file/bslLF165.weights.1 \
    --not-M-5-50 --overlap-annot --print-coefficients \
    --out $PREFIX_OUT

# python /n/groups/price/steven/soft/ldsc/ldsc.py \
#     --h2 $sumstats \
#     --ref-ld /n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001_chimp/ldsc_res/ld_file_c1/bslLF165.1 \
#     --w-ld /n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001_chimp/ldsc_res/ld_file_c1/bslLF165.weights.1 \
#     --not-M-5-50 --overlap-annot --print-coefficients \
#     --out $PREFIX_OUT

# sbatch -p short -t 0-1:00 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $PREFIX_OUT.sbatch.log --wrap " \
# python /n/groups/price/steven/soft/ldsc/ldsc.py \
#     --h2 $sumstats \
#     --ref-ld /n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001_chimp/ldsc_res/ld_file_c1/bslLF165. \
#     --w-ld /n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001_chimp/ldsc_res/ld_file_c1/bslLF165.weights. \
#     --not-M-5-50 --overlap-annot --print-coefficients \
#     --out $PREFIX_OUT"