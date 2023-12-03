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

# Data files
PGEN_FILE=${LDSPEC_PATH}/ldspec/data/data.geno.c@ # `@` refers to CHR number
ANNOT_FILE=${LDSPEC_PATH}/ldspec/data/data.c@.annot.gz
PANNOT_FILE=${LDSPEC_PATH}/ldspec/data/proxy_0_1000.ld_n100_p100.maf_all.c@.pannot_mat.npz
# Output path
OUT_PATH=./results
[ -d ${OUT_PATH} ] || mkdir ${OUT_PATH}


# # Step 1: get SNP blocks
# # get_snp_block : generate SNP blocks (10,000 SNPs per block)
# #     - Input : --job | --pgen_file | --prefix_out
# #     - Output : list of snp ranges (in the format of "snp_range"), e.g., `c1_s0_e10000`.
# python3 ${LDSPEC_PATH}/CLI_ldspec.py\
#     --job get_snp_block\
#     --pgen_file $PGEN_FILE\
#     --prefix_out ${OUT_PATH}/ldspec_s1

# # Step 2: compute LD matrices
# # Compute LD matrix for each SNP block in the snp_range.txt file; can parallelize.
# # compute_ld : compute LD between target SNPs in `snp_range` and reference SNPs within `win_size` of target SNPs
# #     - Input : --job | --pgen_file | --prefix_out | --snp_range | [--win_size]
# #     - Output : LD matrix.
# # '_ld' can only exists in '.<snp_range>_ld.' segment
# SNP_RANGE_FILE=${OUT_PATH}/ldspec_s1.snp_range.txt
# N_SNP_RANGE=$( wc -l < ${SNP_RANGE_FILE} )
# for i_line in $(eval echo "{1..$N_SNP_RANGE}")
# do
# SNP_RANGE=$( head -n $i_line ${SNP_RANGE_FILE} | tail -1 )
# python3 ${LDSPEC_PATH}/CLI_ldspec.py\
#     --job compute_ld\
#     --pgen_file $PGEN_FILE\
#     --snp_range $SNP_RANGE\
#     --win_size 1e6\
#     --prefix_out ${OUT_PATH}/ldspec_s2
# done

# # Step 3: compute LD score
# # Compute LD score for each SNP block in the snp_range.txt file; can parallelize.
# # --annot_file takes comma seperated .annot.gz and .pannot_mat.npz files, or a .txt file of file paths
# SNP_RANGE_FILE=${OUT_PATH}/ldspec_s1.snp_range.txt
# N_SNP_RANGE=$( wc -l < ${SNP_RANGE_FILE} )
# for i_line in $(eval echo "{1..$N_SNP_RANGE}")
# do
# SNP_RANGE=$( head -n $i_line ${SNP_RANGE_FILE} | tail -1 )
# LD_FILE=${OUT_PATH}/ldspec_s2.${SNP_RANGE}_ld.npz
# python3 ${LDSPEC_PATH}/CLI_ldspec.py\
#     --job compute_score\
#     --pgen_file $PGEN_FILE\
#     --ld_file $LD_FILE\
#     --annot_file ${ANNOT_FILE},${PANNOT_FILE}\
#     --win_size 1e6\
#     --prefix_out ${OUT_PATH}/ldspec_s3
# done

# # Step 4: combine LD scores
# #     --score_file: @ represents `snp_range`
# SNP_RANGE_FILE=${OUT_PATH}/ldspec_s1.snp_range.txt
# python3 ${LDSPEC_PATH}/CLI_ldspec.py\
#     --job combine_score\
#     --snp_range_file ${SNP_RANGE_FILE}\
#     --score_file ${OUT_PATH}/ldspec_s3.@_score.tsv.gz\
#     --prefix_out ${OUT_PATH}/ldspec_s4

# # Step 5: compute AVGR
# # compute average LD (avgr) for each pannot. "--ld_file" should contain all LD files
# #     - Input : --job | --pgen_file | --ld_file | --annot_file | --prefix_out
# #     - Output : Average LD for each pannot.
# LD_FILE=${OUT_PATH}/ldspec_s2.@_ld.npz
# python3 ${LDSPEC_PATH}/CLI_ldspec.py\
#     --pgen_file $PGEN_FILE\
#     --job compute_avgr\
#     --ld_file $LD_FILE\
#     --annot_file ${PANNOT_FILE}\
#     --prefix_out ${OUT_PATH}/ldspec_s5

# Step 6: regress
# --avgr_file supports comma-separated .avgr files
SCORE_FILE=${OUT_PATH}/ldspec_s4.c@_score.tsv.gz
AVGR_FILE=${OUT_PATH}/ldspec_s5.avgr
SUMSTATS_FILE=${LDSPEC_PATH}/ldspec/data/simu.seed_0.sumstats.gz
python3 ${LDSPEC_PATH}/CLI_ldspec.py\
    --job regress\
    --pgen_file $PGEN_FILE\
    --annot_file ${ANNOT_FILE},${PANNOT_FILE}\
    --score_file $SCORE_FILE\
    --sumstats_file $SUMSTATS_FILE\
    --avgr_file $AVGR_FILE\
    --prefix_out ${OUT_PATH}/ldspec_s6.seed_0\
    --flag_nofil_snp True    