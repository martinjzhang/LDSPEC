#!/bin/bash

# $ldscore = "python /n/groups/price/steven/soft/ldsc/ldsc.py";
# $options = "--overlap-annot --print-coefficients"; # --print-cov --print-delete-vals";

# $freq    = "--frqfile-chr /n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC.";
# $weights = "--w-ld-chr /n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/weights/weights.hm3_noMHC.";
# $annots1 = "/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baselineLD_v2.2/baselineLD.";

SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.100123/null_h2g50_p20
SIMU_PATH=/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/simulation.100123/causal_h2g50_p20

ANALYSIS=ldsc_debug2
OUT_PATH=${SIMU_PATH}.${ANALYSIS}
[ -d $OUT_PATH ] || mkdir $OUT_PATH

for i_rep in {0..49}
# for i_rep in 1
do
sumstats=${SIMU_PATH}/rep${i_rep}.sumstats.gz
PREFIX_OUT=${OUT_PATH}/rep${i_rep}

sbatch -p short -t 0-0:20 -n 1 -c 1 --mem=16000 --open-mode=truncate -o $PREFIX_OUT.sbatch.log --wrap " \
python /n/groups/price/steven/soft/ldsc/ldsc.py \
    --h2 $sumstats \
    --ref-ld /n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldsc_res/ld_file/bslLF165.1 \
    --w-ld /n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001_chimp/ldsc_res/ld_file/bslLF165.weights.1 \
    --not-M-5-50 --overlap-annot --print-coefficients \
    --out $PREFIX_OUT"
done

#     --chisq-max 100000000 \
