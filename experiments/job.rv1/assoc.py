import pandas as pd
import numpy as np
import scipy as sp
import os
from os.path import join
import re
import time
import argparse
import pickle
import ldspec


"""
Job description
----------------

get_snp_block : create a list of SNP blocks (10,000 SNPs per block)
    - Input : --job | --pgen_file | --prefix_out
    - Output : one line per block in the `snp_range` format, e.g., `c1_s0_e10000`.

compute_ld : compute LD between target SNPs in `snp_range` and reference SNPs within `win_size` of target SNPs
    - Input : --job | --pgen_file | --prefix_out | --snp_range | [--win_size]
    - Output : `_ld.npz` file; LD matrix.
    
compute_score : compute LD and DLD scores.
    - Input : --job | --pgen_file | --ld_file | --annot_file | --prefix_out | [--win_size] | [--flag_cross_term]
    - Output : `_score.tsv.gz` file; LD and DLD scores.
    
combine_score : concatenate scores from the same CHR 
    - Input : --score_file | --snp_range_file | --prefix_out
    - Output : concatenated score files by CHR
    
compute_avgr : compute average LD (avgr) for each pannot. "--ld_file" should contain all LD files
    - Input : --job | --pgen_file | --ld_file | --annot_file | --prefix_out
    - Output : Average LD for each pannot.
    
regress : estimate LDSPEC parameters.
    - Input : --job | --pgen_file | --annot_file  | --score_file | --sumstats_file| --avgr_file | --prefix_out |
    [--flag_cross_term] | [--flag_nofil_snp]
    - Output : LD-SPEC result.
    
evaluate : model evaluation 
    - Input : --job | --pgen_file | --annot_file  | --score_file | --sumstats_file| --avgr_file | --null_model_file | 
    --prefix_out | [--flag_cross_term] | [--flag_nofil_snp]
    - Output : model evaluation results.
    
TODO
----
- ldspec.util.update_columns ALT_FREQS is updated as MAF, which is not correct for derived alleles. It may not cause any trouble now as MAF is used in a sysmetrical way: p (1-p). But a thorough check is needed.
- compute_score: output SNP alignment results
- 
"""


def main(args):
    sys_start_time = time.time()

    ###########################################################################################
    ######                                    Parse Options                              ######
    ###########################################################################################

    TRAIT = args.trait
    PREFIX_OUT = args.prefix_out
    
    # Print input options
    header = ldspec.util.get_cli_head()
    header += "Call: assoc.py \\\n"
    header += "--trait %s\\\n" % TRAIT
    header += "--prefix_out %s\\\n" % PREFIX_OUT
    print(header)
    
    trait = TRAIT
    CHR_LIST = np.arange(1,23).tolist()
    # CHR_LIST = [1, 7, 10, 21]
    ###########################################################################################
    ######                                  Data Loading                                 ######
    ###########################################################################################
    DATA_PATH = "/n/groups/price/martin/LDSPEC_data/"
    sumstats_file = '/n/groups/price/martin/LDSPEC_data/UKBimp_ss50K_MAF001_chimp/sumstats/@.sumstats.gz'
    trait_list = pd.read_csv(
        '/n/groups/price/martin/LDSPEC_data/UKBB_trait_ss50K/trait_list_all_indpt.txt', header=None,
    )[0].to_list()

    # genotype files
    pgen_file = '/n/scratch/users/j/jz286/imp_geno_ss50K_chimp/ukb_imp_chr@_v3_chimp'
    dic_data = {}
    for CHR in range(1, 23):  # Check all 23 CHRs
        if os.path.exists(pgen_file.replace("@", "%s" % CHR) + ".pgen"):
            dic_data[CHR] = ldspec.util.read_pgen(
                pgen_file.replace("@", "%s" % CHR)
            )
    ld_file = '/n/scratch/users/j/jz286/imp_geno_ss50K_chimp.ldspec_ld_1e6/ukb_imp_v3_chimp.@_ld.npz'
    snp_range_list = pd.read_csv(
        '/n/groups/price/martin/LDSPEC_data/UKBimp_ss50K_MAF001_chimp/ukb_imp_v3.snp_range.txt', header=None,
    )[0].to_list()
    
    df_sumstats = pd.read_csv(sumstats_file.replace('@', trait), sep="\t", index_col=None)
    df_sumstats.index = df_sumstats['SNP']
    
    # df_snp
    df_snp = []
    for CHR in CHR_LIST:
        df_snp_chr = dic_data[CHR]['pvar'].copy()
        df_snp_chr.index = df_snp_chr['SNP']
        temp_dic = {x:y for x,y in zip(dic_data[CHR]['afreq']['SNP'],dic_data[CHR]['afreq']['MAF'])}
        df_snp_chr['AF'] = [temp_dic[x] if x in temp_dic else 0 for x in df_snp_chr['SNP']]
        # Add SNP range    
        df_snp_chr['idx'] = np.arange(df_snp_chr.shape[0])
        df_snp_chr['idx_range'] = np.arange(df_snp_chr.shape[0])
        df_snp_chr['range'] = ''
        for temp_range in snp_range_list:
            if temp_range.startswith('c%d_'%CHR):
                temp_dic = ldspec.util.parse_snp_range(temp_range)
                ind_select = np.zeros(df_snp_chr.shape[0], dtype=bool)
                ind_select[temp_dic['start']:temp_dic['end']] = True
                df_snp_chr.loc[ind_select, 'range'] = temp_range
                df_snp_chr.loc[ind_select, 'idx_range'] -= temp_dic['start']
        df_snp.append(df_snp_chr)
    df_snp = pd.concat(df_snp, axis=0)
    # add sumstats 
    temp_dic = {x:y for x,y in zip(df_sumstats['SNP'],df_sumstats['Z'])}
    df_snp['Z'] = [temp_dic[x] if x in temp_dic else 0 for x in df_snp['SNP']]
    df_snp['P'] = ldspec.util.zsc2pval(df_snp['Z'])
    print("    df_snp loaded ", df_snp.shape)
    print("    " + ldspec.util.get_sys_info(sys_start_time))

    ###########################################################################################
    ######                                  Computation                                  ######
    ###########################################################################################
    # Pairwise testing: P<5e-6 SNPs and other SNPs within 10kb and P<1e-3
    # Pairwise testing: P<5e-4 SNPs and other SNPs within 10kb and P<5e-2
    # (snp1 has a smaller p-value than snp2 to avoid duplicates)
#     snp_analyze_list = df_snp['SNP'][df_snp['P']<5e-6].to_list()
    snp_analyze_list = df_snp['SNP'][df_snp['P']<5e-4].to_list()
    # Testing 
    rho_r_ratio_list = [-1, -0.5, 0, 0.5, 1]
    win_size = 1e4
    dic_res = {x:[] for x in rho_r_ratio_list}

    # snp_analyze_list = snp_analyze_list[:100]
    # rho_r_ratio_list = [-1, -0.5, 0]

    df_snp_analyze = df_snp.loc[snp_analyze_list] # info for the set of SNPs to analyze
    print('Analyzing %d SNPs in %d range files' % (
        df_snp_analyze.shape[0], df_snp_analyze['range'].unique().shape[0]
    ))
    for temp_range in df_snp_analyze['range'].unique():
        # One LD file at a time
        print('    ', temp_range)
        temp_dic = ldspec.util.parse_snp_range(temp_range)
        start,end = temp_dic['start'],temp_dic['end']
        mat_ld, dic_range = ldspec.util.read_ld(ld_file.replace('@', temp_range))
        df_snp_chr = df_snp.loc[df_snp['CHR']==temp_dic['chr']] # SNPs on the same CHR
        input_list = []

        for snp1 in df_snp_analyze.loc[df_snp_analyze['range']==temp_range, 'SNP']:
            temp_bp = df_snp_chr.loc[snp1, 'BP']
            v_bp = df_snp_chr['BP']
            v_p = df_snp_chr['P']
            # within 10kb
            # ind_select = (v_bp>temp_bp-win_size) & (v_bp<temp_bp+win_size) & (v_p>df_snp_chr.loc[snp1, 'P']) & (v_p<1e-3)
            ind_select = (v_bp>temp_bp-win_size) & (v_bp<temp_bp+win_size) & (v_p>df_snp_chr.loc[snp1, 'P']) & (v_p<5e-2)
            for snp2 in df_snp_chr.loc[ind_select, 'SNP']:
                ld = mat_ld[df_snp_chr.loc[snp2, 'idx'], df_snp_chr.loc[snp1, 'idx_range']]
                abs_ld = np.absolute(ld)
                #  abs(LD)>0.1, abs(LD)<0.9, P>P_focal
                #  abs(LD)>0.05, abs(LD)<0.95, P>P_focal
                if (abs_ld>0.05) & (abs_ld<0.95):
                    input_list.append([snp1, snp2, df_snp_chr.loc[snp1, 'Z'], df_snp_chr.loc[snp2, 'Z'], ld, 5e4]) 
        for rho_r_ratio in rho_r_ratio_list:
            temp_df = ldspec.util.analyze_batch(input_list, rho_r_ratio=rho_r_ratio)
            dic_res[rho_r_ratio].append(temp_df)

    for rho_r_ratio in rho_r_ratio_list:
        dic_res[rho_r_ratio] = pd.concat(dic_res[rho_r_ratio], axis=0)
        dic_res[rho_r_ratio].to_csv(PREFIX_OUT+'.%s.%.1f.tsv' % (trait, rho_r_ratio), sep="\t", index=False)
    print("    " + ldspec.util.get_sys_info(sys_start_time))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ldspec")
    parser.add_argument("--trait", type=str, required=True, help="")
    parser.add_argument("--prefix_out", type=str, required=True)
    args = parser.parse_args()
    main(args)
