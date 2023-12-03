#!/bin/bash

# $ldscore = "python /n/groups/price/steven/soft/ldsc/ldsc.py";
# $options = "--overlap-annot --print-coefficients"; # --print-cov --print-delete-vals";

# $freq    = "--frqfile-chr /n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC.";
# $weights = "--w-ld-chr /n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/weights/weights.hm3_noMHC.";
# $annots1 = "/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baselineLD_v2.2/baselineLD.";


for i_line in {2..29}
# for i_line in 1
do
TRAIT=$( head -n $i_line "/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001/UKBB_trait/trait_list_all.txt" | tail -1 )
# TRAIT=blood_PLATELET_COUNT
sumstats=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldsc_res/sumstats/${TRAIT}.nomhc.sumstats.gz
PATH_OUT=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldsc_res
PREFIX_OUT=${PATH_OUT}/${TRAIT}.bslLF165

sbatch -p short -t 0-1:00 -n 1 -c 1 --mem=128000 --open-mode=truncate -o $PREFIX_OUT.sbatch.log --wrap " \
python /n/groups/price/steven/soft/ldsc/ldsc.py \
    --h2 $sumstats \
    --ref-ld-chr /n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldsc_res/ld_file/bslLF165. \
    --w-ld-chr /n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldsc_res/ld_file/bslLF165.weights. \
    --not-M-5-50 --overlap-annot --print-coefficients \
    --out $PREFIX_OUT"
done