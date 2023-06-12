#SBATCH -c 8
#SBATCH -N 1
#SBATCH -t 0-01:00
#SBATCH --array=1-7
#SBATCH -p short
#SBATCH --mem=64000
#SBATCH -o /home/jz286/WES_analysis/wes_rare/experiments/job_info/hostname_%j.out  
#SBATCH -e /home/jz286/WES_analysis/wes_rare/experiments/job_info/hostname_%j.err 
#SBATCH --mail-type=NONE#SBATCH --mail-type=NONE

# Your local LDSPEC path
LDSPEC_PATH=/home/jz286/WES_analysis/LDSPEC
PGEN_FILE=${LDSPEC_PATH}/ldspec/data/data.geno.c@ # `@` refers to CHR number
ANNOT_FILE=${LDSPEC_PATH}/ldspec/data/data.c@.annot.gz
PANNOT_FILE=${LDSPEC_PATH}/ldspec/data/proxy_0_1000.ld_n100_p100.maf_all.c@.pannot_mat.npz

# Simulate data: 
#     SNP effect file `.eff.gz`
#     Phenotype file `.phen`
#     Summary statistics file `.sumstats.gz`
#     Effect summary file `.eff_tau.tsv` and `.eff_omega.tsv`
# Example config file: (.tsv file, `h2g`, `p_causal`, `alpha` required)
#     h2g     0.500
#     p_causal        0.900
#     alpha   -1
#     AN:all  1.000
#     AN:DHS_Trynka_common        1.000
#     pAN:proxy_0_1000_ld_n100_p100_maf_all   0.500
# `alpha`: Parameter for maf-dependent architecture, `Var(beta_j) \propto [maf (1-maf)]^(1+\alpha)`, 
# where `\beta_j` is the standardized SNP effect size.  `alpha=-1` means no MAF-dependency. 
# Schoech NC 2019 suggested alpha=-0.38.     
CONFIG_FILE=${LDSPEC_PATH}/ldspec/data/config.tsv
SEED=0
python3 ${LDSPEC_PATH}/CLI_simulation.py\
    --job simulate\
    --pgen_file $PGEN_FILE\
    --config_file $CONFIG_FILE\
    --annot_file ${ANNOT_FILE},${PANNOT_FILE}\
    --random_seed ${SEED}\
    --prefix_out ${LDSPEC_PATH}/ldspec/data/simu.seed_${SEED}