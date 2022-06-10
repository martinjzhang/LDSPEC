import numpy as np
import scipy as sp
import pandas as pd
import time
import gdreg
import warnings


def estimate(
    df_score
    df_sumstats,
    dic_data,
    df_annot,
    pannot_list=[],
    pannot_hr_list=[],
    n_jn_block=100,
    sym_non_pAN="non-pAN",
    win_size=int(1e7),
    memory=128,
    n_iter=1,
    verbose=False,
):
    """
    Multi-stage GDREG estimation.

    Parameters
    ----------
    df_score : pd.DataFrame, default=None
        GDReg LD and DLD scores, with columns ['CHR', 'SNP', 'BP', 'LD:AN1',
        'LD:AN2', 'DLD:PAN:AN1', 'DLD:PAN:AN2'].
    dic_data : dict
        Genotype data reader, organized, for each CHR, as

        - dic_data[CHR]['pgen'] : .pgen file path
        - dic_data[CHR]['pvar'] : .pvar pd.DataFrame
        - dic_data[CHR]['psam'] : .psam pd.DataFrame

    df_sumstats : pd.DataFrame
        Summary statistics with columns ['SNP', 'N', 'Z', 'A1', 'A2']
    dic_ld : dict of sp.sparse.csr_matrix (TODO)
        dic_mat_ld[CHR] for LD matrix of chromosome CHR. Used in `gdreg.score.compute_score`.
    df_annot : pd.DataFrame
        Single-SNP annotation. Used in `gdreg.score.summarize`.
    df_pannot_list : list of pd.DataFrame, default=[]
        Each element corresponds to SNP-pair annotation. Must contain
        ['CHR', 'SNP', 'BP', 'pAN:pAN1'] columns. Used in `gdreg.score.compute_score`.
    df_pannot_hr_list : list of pd.DataFrame, default=[]
        Each element corresponds to a high-res SNP-pair annotation. Must contain
        ['CHR', 'SNP', 'BP', 'pCHR', 'pSNP', 'pBP', 'pAN:pAN1'] columns. Used in
        `gdreg.score.compute_score`.
    n_jn_block : int
        Number of JN blocks.
    sym_non_pAN : str, default='non-pAN'
        Symbol for SNPs not in the SNP-pair annotation.
    win_size : int, defualt=1e7
        Window size for computing LD and DLD scores. Used in `gdreg.score.compute_score`.
    memory : int, default=128
        Memory to use (in MB). Used in `gdreg.score.compute_score`.
    verbose : bool, default=False
        If to output messages.

    Returns
    -------
    df_res : pd.DataFrame
        GDREG regression results.
    """

    start_time = time.time()
    if verbose:
        print("# Call: gdreg.regress.estimate")

    # df_snp from dic_data
    CHR_list = sorted(dic_data)
    df_snp = None
    for CHR in CHR_list:
        if df_snp is None:
            df_snp = dic_data[CHR]["pvar"].copy()
        else:
            df_snp = pd.concat([df_snp, dic_data[CHR]["pvar"]], axis=0)

    if verbose:
        print(
            "    dic_data : n_snp=%d, n_sample=%d" % (df_snp.shape[0], df_snp.shape[0])
        )

    # df_sumstats
    n_sample_zsq = df_sumstats["N"].mean().astype(int)
    dic_zsq = {x: y**2 for x, y in zip(df_sumstats["SNP"], df_sumstats["Z"])}
    outlier_thres = max(80, 0.001 * n_sample_zsq)  # Finucane 2015 Nat Genet
    dic_zsq = {x: y for x, y in dic_zsq.items() if y < outlier_thres}

    if verbose:
        print(
            "    df_sumstats : n_snp=%d, n_sample_zsq=%d"
            % (df_sumstats.shape[0], n_sample_zsq)
        )
        print(
            "    Remove duplicates or ZSQ>%0.1f SNPs, %d remaining, avg. zsq=%0.2f"
            % (outlier_thres, len(dic_zsq), np.mean(list(dic_zsq.values())))
        )

    # df_reg
    df_reg = df_snp.copy()
    df_reg.drop_duplicates("SNP", inplace=True)
    df_reg = df_reg.loc[df_reg["SNP"].isin(dic_zsq)]
    df_reg["ZSQ"] = [dic_zsq[x] for x in df_reg["SNP"]]
    df_reg.index = df_reg["SNP"]
    df_reg.sort_values(by=["CHR", "BP"], inplace=True)
    dic_block = get_block(df_reg, pannot_list, sym_non_pAN, n_jn_block)

    if verbose:
        print(
            "    Regression : n_snp=%d, n_block=%d" % (df_reg.shape[0], len(dic_block))
        )

    # Regression : LD-score only and estimate \tau
    dic_res = {}

    if df_score is None:
        df_score = gdreg.score.compute_score(
            dic_data,
            dic_ld,
            df_annot,
            pannot_list=pannot_list,
            pannot_hr_list=pannot_hr_list,
            sym_non_pAN=sym_non_pAN,
            win_size=win_size,
            memory=memory,
            verbose=verbose,
            verbose_prefix='    ',
        )
#         df_score = gdreg.score.compute_score(
#             dic_data,
#             dic_ld,
#             df_annot,
#             win_size=win_size,
#             memory=memory,
#             verbose=verbose,
#         )

    df_score.index = df_score["SNP"]
    reg_list = [x for x in df_score if x.startswith("LD:")] + ["E"]
    temp_df_reg = df_reg.join(df_score[reg_list])
    dic_res[0] = regress(temp_df_reg, dic_block, n_sample_zsq, verbose=verbose, verbose_prefix='    ')
#     dic_res[0]["summary"] = summarize(dic_res[0], df_annot)
    
    # Regression : estimate both \tau and \rho
    reg_list = [x for x in df_score if x.startswith("LD:")] + [x for x in df_score if x.startswith("DLD:")] + ["E"]
    temp_df_reg = df_reg.join(df_score[reg_list])
    dic_res[1] = regress(temp_df_reg, dic_block, n_sample_zsq, verbose=verbose, verbose_prefix='    ')
#     dic_res[1]["summary"] = summarize(dic_res[1], df_annot)
    
#     # Second pass : \tau & \rho
#     for i_iter in range(1, n_iter+1):
#         v_var_ps = np.zeros(df_annot.shape[0], dtype=np.float32)
#         res = i_iter - 1
#         for LDAN, coef in zip(dic_res[res]["term"], dic_res[res]["coef_jn"]):
#             if (LDAN != "E") & (LDAN.startswith("DLD:") is False):
#                 AN = LDAN.replace("LD:", "")
#                 v_var_ps += coef * df_annot[AN].values
#         v_var_ps = v_var_ps.clip(
#             min=0.001 / v_var_ps.shape[0]
#         )  # Lower bound (check later)
#         dic_var_ps = {x: y for x, y in zip(df_annot["SNP"], v_var_ps)}

#         # Recompute score (computational bottleneck)
#         df_score = gdreg.score.compute_score(
#             dic_data,
#             dic_ld,
#             df_annot,
#             pannot_list=pannot_list,
#             pannot_hr_list=pannot_hr_list,
#             sym_non_pAN=sym_non_pAN,
#             win_size=win_size,
#             memory=memory,
#             verbose=verbose,
#             verbose_prefix='    ',
#             dic_var_ps=dic_var_ps,
#         )

#         df_score.index = df_score["SNP"]
#         reg_list = (
#             [x for x in df_score if x.startswith("LD:")]
#             + [x for x in df_score if x.startswith("DLD:")]
#             + ["E"]
#         )
#         temp_df_reg = df_reg.join(df_score[reg_list])
#         dic_res[i_iter] = regress(
#             temp_df_reg, dic_block, n_sample_zsq, verbose=verbose, verbose_prefix='    ',
#         )
        
#         dic_res[i_iter]["summary"] = summarize(dic_res[i_iter], df_annot)
        
#         # Convergence check
#         v_tau = dic_res[i_iter]["summary"]["tau"]
#         v_tau_prev = dic_res[i_iter-1]["summary"]["tau"]
#         rel_dif = np.absolute(v_tau-v_tau_prev).sum() / np.absolute(v_tau).sum()
#         dic_res[i_iter]["converge"] = rel_dif<0.01 
#         if verbose:
#             print("    Iter %d finished, ref_dif=%0.2g, converge=%s\n" % (i_iter, rel_dif, dic_res[i_iter]["converge"]))
#         if dic_res[i_iter]["converge"]:
#             break
            
    if verbose:
        print("    Completed, time=%0.1fs" % (time.time() - start_time))

    return dic_res


def summarize(dic_res, df_annot,
              pannot_list=[],
              pannot_hr_list=[],
              sym_non_pAN="non-pAN",
              ):
    """
    Summarize GDREG result.

    Parameters
    ----------
    dic_res : dict
        Regression results.

        - dic_res['term'] : list of terms.
        - dic_res['coef'] : estimated coefs. np.ndarray(dtype=np.float32).
        - dic_res['coef_jn'] : JN estimated coefs. np.ndarray(dtype=np.float32).
        - dic_res['coef_jn_cov'] : estimated coef covariance. np.ndarray(dtype=np.float32).

    df_annot : pd.DataFrame
        Single-SNP annotation. Used in `gdreg.score.compute_score`.
    df_pannot_list : list of pd.DataFrame, default=[]
        Each element corresponds to SNP-pair annotation. Must contain
        ['CHR', 'SNP', 'BP', 'pAN:pAN1'] columns. Used in `gdreg.score.compute_score`.
    df_pannot_hr_list : list of pd.DataFrame, default=[]
        Each element corresponds to a high-res SNP-pair annotation. Must contain
        ['CHR', 'SNP', 'BP', 'pCHR', 'pSNP', 'pBP', 'pAN:pAN1'] columns. Used in
        `gdreg.score.compute_score`.
    sym_non_pAN : str, default='non-pAN'
        Symbol for SNPs not in the SNP-pair annotation.
    n_block: int, default=100
        Number of jackknife blocks.

    Returns
    -------
    df_summary : pd.DataFrame
        Regression result summary.
        
    TODO
    ----
    1. Double check later.

    """
    
    # Annotation info
    AN_list = [x for x in df_annot if x.startswith("AN:")]
    pAN_list = [[y for y in x if y.startswith("pAN")][0] for x in pannot_list]
    pAN_hr_list = [[y for y in x if y.startswith("pAN")][0] for x in pannot_hr_list]

    dic_annot = {}
    for AN in AN_list:
        dic_annot[AN] = {x: y for x, y in zip(df_annot["SNP"], df_annot[AN])}
    dic_pannot = {}
    for temp_df in pannot_list:
        pAN = [x for x in temp_df if x.startswith("pAN:")][0]
        dic_pannot[pAN] = {x: y for x, y in zip(temp_df["SNP"], temp_df[pAN])}
    dic_pannot_hr = {}
    for temp_df in pannot_hr_list:
        pAN = [x for x in temp_df if x.startswith("pAN:")][0]
        dic_pannot_hr[pAN] = [(x, y) for x, y in zip(temp_df["SNP"], temp_df["pSNP"])]
        
    # dic_mat_G
    dic_mat_G = {}
    v_snp = np.array(df_annot["SNP"])
    snp_set = set(v_snp)
    for pAN in pAN_list:       
        v_pAN = [
                dic_pannot[pAN][x] if x in dic_pannot[pAN] else sym_non_pAN
                for x in v_snp
            ]
        mat_G = gdreg.util.pannot_to_csr(v_pAN, sym_non_pAN=sym_non_pAN, flag_matS=False)
        dic_mat_G[pAN] = mat_G.copy()

    for pAN in pAN_hr_list:
        snp_pair_list = [
            x
            for x in dic_pannot_hr[pAN]
            if (x[0] in snp_set) & (x[1] in snp_set)
        ]
        mat_G = gdreg.util.pannot_hr_to_csr(v_snp, snp_pair_list)
        dic_mat_G[pAN] = mat_G.copy()    

    # Results info
    LD_list = [x for x in dic_res["term"] if x.startswith("LD:")]
    DLD_list = [x for x in dic_res["term"] if x.startswith("DLD:")]

    res_AN_list = [x.replace("LD:", "") for x in LD_list]
    res_pAN_list = sorted(set([x.replace("DLD:", "").split("|")[0] for x in DLD_list]))

    # Check consistency betweeno results and annotation df_annot
    err_msg = "df_annot does not contain all annots in dic_res"
    assert len(set(res_AN_list) - set(AN_list)) == 0, err_msg 
    err_msg = "pannot_list and pannot_hr_list does not contain all pannots in dic_res"
    assert len(set(res_pAN_list) - set(pAN_list+pAN_hr_list)) == 0, err_msg 
    
    col_list = ["annot", "n_snp", "tau", "tau_se", "h2", "h2_se", "enrich", "enrich_se"]
    df_summary = pd.DataFrame(index=res_AN_list, columns=col_list)
    col_list = ["pannot", "n_snp_pair", "rho", "rho_se", "cov", "cov_se"]
    df_summary_p = pd.DataFrame(index=res_pAN_list, columns=col_list)
    
    # Summary : basic
    df_summary["annot"] = df_summary.index
    df_summary["n_snp"] = [(df_annot[x] == 1).sum() for x in res_AN_list]
    df_summary_p["pannot"] = df_summary_p.index
    df_summary_p["n_snp_pair"] = [dic_mat_G[x].sum() for x in res_pAN_list]

    # Summary : tau & rho
    temp_dic = {x: y for x, y in zip(dic_res["term"], dic_res["coef_jn"])}
    df_summary["tau"] = [temp_dic["LD:%s" % x] for x in res_AN_list]
    df_summary_p["rho"] = [temp_dic["DLD:%s" % x] for x in res_pAN_list]

    temp_dic = {x: np.sqrt(y) for x, y in zip(dic_res["term"], np.diag(dic_res["coef_jn_cov"]))}
    df_summary["tau_se"] = [temp_dic["LD:%s" % x] for x in res_AN_list]
    df_summary_p["rho_se"] = [temp_dic["DLD:%s" % x] for x in res_pAN_list]

    # Summary : h2 & h2_se & enrich
    df_cov = pd.DataFrame(
        index=dic_res["term"], columns=dic_res["term"], data=dic_res["coef_jn_cov"]
    )
    temp_mat = df_cov.loc[LD_list, LD_list].values
    for AN in res_AN_list:
        if len(set(df_annot[AN])) > 2:
            continue
            
        temp_v = df_annot.loc[df_annot[AN] == 1, res_AN_list].sum(axis=0).values
        df_summary.loc[AN, "h2"] = (temp_v * df_summary["tau"]).sum()
        df_summary.loc[AN, "h2_se"] = np.sqrt(temp_v.dot(temp_mat).dot(temp_v))

        # enrich
        # determine reference annotation
        # TODO : test
        ind_ref = np.ones(df_annot.shape[0], dtype=bool)
        for term in ["_common", "_lf", "_rare"]:
            if AN.endswith(term):
                col_list = [x for x in df_annot if x.endswith(term)]
                ind_ref = (df_annot[col_list].values == 1).sum(axis=1) > 0
        temp_v_ref = df_annot.loc[ind_ref, res_AN_list].sum(axis=0).values

        n_snp_cate = (df_annot[AN] == 1).sum()
        n_snp_dif = ind_ref.sum() - n_snp_cate
        if n_snp_dif < n_snp_cate*0.1:
            continue
        temp_v_combine = (
            temp_v * (1 / n_snp_cate + 1 / n_snp_dif) - temp_v_ref / n_snp_dif
        )

        h2_ps_ref = (temp_v_ref * df_summary["tau"]).sum() / ind_ref.sum()
        df_summary.loc[AN, "enrich"] = df_summary.loc[AN, "h2"] / n_snp_cate / h2_ps_ref

        dif_ = (temp_v_combine * df_summary["tau"]).sum()
        se_ = np.sqrt(temp_v_combine.dot(temp_mat).dot(temp_v_combine))
        df_summary.loc[AN, "enrich_se"] = np.absolute(
            se_ / dif_ * (df_summary.loc[AN, "enrich"] - 1)
        )

    # Summary : r2 and per-pannot r2
    for pAN in res_pAN_list:
        
        temp_v = np.array([(dic_mat_G[pAN] * dic_mat_G[x]).sum() for x in res_pAN_list])
        temp_rho = df_summary_p['rho'].values
        temp_mat = df_cov.loc[DLD_list, DLD_list].values
        df_summary_p.loc[pAN, "cov"] = (temp_v * temp_rho).sum()
        df_summary_p.loc[pAN, "cov_se"] = np.sqrt(temp_v.dot(temp_mat).dot(temp_v))

    return df_summary,df_summary_p


# def summarize(dic_res, df_annot):
#     """
#     Summarize GDREG result.

#     Parameters
#     ----------
#     dic_res_reg : dict
#         Regression results.

#         - dic_res_reg['term'] : list of terms.
#         - dic_res_reg['coef'] : estimated coefs. np.ndarray(dtype=np.float32).
#         - dic_res_reg['coef_jn'] : JN estimated coefs. np.ndarray(dtype=np.float32).
#         - dic_res_reg['coef_jn_cov'] : estimated coef covariance. np.ndarray(dtype=np.float32).

#     df_annot : pd.DataFrame
#         Single-SNP annotation. Used in `gdreg.score.compute_score`.
#     df_pannot_list : list of pd.DataFrame, default=[]
#         Each element corresponds to SNP-pair annotation. Must contain
#         ['CHR', 'SNP', 'BP', 'pAN:pAN1'] columns. Used in `gdreg.score.compute_score`.
#     df_pannot_hr_list : list of pd.DataFrame, default=[]
#         Each element corresponds to a high-res SNP-pair annotation. Must contain
#         ['CHR', 'SNP', 'BP', 'pCHR', 'pSNP', 'pBP', 'pAN:pAN1'] columns. Used in
#         `gdreg.score.compute_score`.
#     sym_non_pAN : str, default='non-pAN'
#         Symbol for SNPs not in the SNP-pair annotation.
#     n_block: int, default=100
#         Number of jackknife blocks.

#     Returns
#     -------
#     df_summary : pd.DataFrame
#         Regression result summary.
        
#     TODO
#     ----
#     1. Double check later.

#     """

#     # Annotation info
#     LD_list = [x for x in dic_res["term"] if x.startswith("LD:")]
#     DLD_list = [x for x in dic_res["term"] if x.startswith("DLD:")]

#     AN_list = [x.replace("LD:", "") for x in LD_list]
#     pAN_list = sorted(set([x.replace("DLD:", "").split("|")[0] for x in DLD_list]))

#     # Check consistency with df_annot
#     err_msg = "df_annot does not contain all annots in dic_res"
#     assert len(set(AN_list) - set(df_annot)) == 0, err_msg 

#     col_list = ["annot", "n_snp", "tau", "tau_se", "h2", "h2_se", "enrich", "enrich_se", "r2", "r2_se"]
#     for pAN in pAN_list:
#         col_list += ["%s|%s" % (x, pAN) for x in ["rho", "rho_se", "r2", "r2_se"]]

#     df_summary = pd.DataFrame(index=AN_list, columns=col_list)

#     # Summary : basic
#     df_summary["annot"] = df_summary.index
#     df_summary["n_snp"] = [(df_annot[x] == 1).sum() for x in AN_list]

#     # Summary : tau & rho
#     temp_dic = {x: y for x, y in zip(dic_res["term"], dic_res["coef_jn"])}
#     df_summary["tau"] = [temp_dic["LD:%s" % x] for x in AN_list]
#     for pAN in pAN_list:
#         df_summary["rho|%s" % pAN] = [temp_dic["DLD:%s|%s" % (pAN, x)] for x in AN_list]

#     temp_dic = {x: np.sqrt(y) for x, y in zip(dic_res["term"], np.diag(dic_res["coef_jn_cov"]))}
#     df_summary["tau_se"] = [temp_dic["LD:%s" % x] for x in AN_list]
#     for pAN in pAN_list:
#         df_summary["rho_se|%s" % pAN] = [temp_dic["DLD:%s|%s" % (pAN, x)] for x in AN_list]

#     # Summary : h2 & h2_se & enrich
#     df_cov = pd.DataFrame(
#         index=dic_res["term"], columns=dic_res["term"], data=dic_res["coef_jn_cov"]
#     )
#     temp_mat = df_cov.loc[LD_list, LD_list].values
#     for AN in AN_list:
#         if len(set(df_annot[AN])) > 2:
#             continue
            
#         temp_v = df_annot.loc[df_annot[AN] == 1, AN_list].sum(axis=0).values
#         df_summary.loc[AN, "h2"] = (temp_v * df_summary["tau"]).sum()
#         df_summary.loc[AN, "h2_se"] = np.sqrt(temp_v.dot(temp_mat).dot(temp_v))

#         # enrich
#         # determine reference annotation
#         # TODO : test
#         ind_ref = np.ones(df_annot.shape[0], dtype=bool)
#         for term in ["_common", "_lf", "_rare"]:
#             if AN.endswith(term):
#                 col_list = [x for x in df_annot if x.endswith(term)]
#                 ind_ref = (df_annot[col_list].values == 1).sum(axis=1) > 0
#         temp_v_ref = df_annot.loc[ind_ref, AN_list].sum(axis=0).values

#         n_snp_cate = (df_annot[AN] == 1).sum()
#         n_snp_dif = ind_ref.sum() - n_snp_cate
#         if n_snp_dif < n_snp_cate*0.1:
#             continue
#         temp_v_combine = (
#             temp_v * (1 / n_snp_cate + 1 / n_snp_dif) - temp_v_ref / n_snp_dif
#         )

#         h2_ps_ref = (temp_v_ref * df_summary["tau"]).sum() / ind_ref.sum()
#         df_summary.loc[AN, "enrich"] = df_summary.loc[AN, "h2"] / n_snp_cate / h2_ps_ref

#         dif_ = (temp_v_combine * df_summary["tau"]).sum()
#         se_ = np.sqrt(temp_v_combine.dot(temp_mat).dot(temp_v_combine))
#         df_summary.loc[AN, "enrich_se"] = np.absolute(
#             se_ / dif_ * (df_summary.loc[AN, "enrich"] - 1)
#         )

#     # Summary : r2 and per-pannot r2
#     for AN in AN_list:
#         if len(set(df_annot[AN])) > 2:
#             continue
#         temp_df = df_annot.loc[df_annot[AN] == 1]
#         # r2 
#         pAN_AN_list = [x.replace('DLD:','').split('|') for x in DLD_list]
#         temp_v = np.array([temp_df[x[1]].mean() for x in pAN_AN_list])
#         temp_rho = np.array([df_summary.loc[x[1], "rho|%s"%x[0]] for x in pAN_AN_list])
#         temp_mat = df_cov.loc[DLD_list, DLD_list].values
#         df_summary.loc[AN, "r2"] = (temp_v * temp_rho).sum()
#         df_summary.loc[AN, "r2_se"] = np.sqrt(temp_v.dot(temp_mat).dot(temp_v))
        
#         # per-pannot r2
#         for pAN in pAN_list:
#             sub_DLD_list = [x for x in DLD_list if pAN in x]
#             pAN_AN_list = [x.replace('DLD:','').split('|') for x in sub_DLD_list]
#             temp_v = np.array([temp_df[x[1]].mean() for x in pAN_AN_list])
#             temp_rho = np.array([df_summary.loc[x[1], "rho|%s"%x[0]] for x in pAN_AN_list])
#             temp_mat = df_cov.loc[sub_DLD_list, sub_DLD_list].values
#             df_summary.loc[AN, "r2|%s"%pAN] = (temp_v * temp_rho).sum()
#             df_summary.loc[AN, "r2_se|%s"%pAN] = np.sqrt(temp_v.dot(temp_mat).dot(temp_v))

#     return df_summary


# https://github.com/martinjzhang/wes_rare/blob/master/wes_rare/ld_method.py


def get_block(df_reg, pannot_list=[], sym_non_pAN="non-pAN", n_block=100):
    """
    Block partition. SNPs on the same .pannot annotation go to the same block.
    SNPs on different CHRs go to different blocks.

    Parameters
    ----------
    df_reg : pd.DataFrame
        SNPs used in regression. Assumed to be sorted by genomic location.
        Should contain ['CHR', 'SNP', 'BP', 'ZSQ'].
    pannot_list : list of pd.DataFrame, default=[]
        Each element corresponds to SNP-pair annotation. Must contain
        ['CHR', 'SNP', 'BP', 'pAN:pAN1'] columns.
    sym_non_pAN : str, default='non-pAN'
        Symbol for SNPs not in the SNP-pair annotation.
    n_block: int, default=100
        Number of jackknife blocks.

    Returns
    -------
    dic_block : dict
        Block information. `dic_block[i] = (ind_s, ind_e)`.

    """

    n_snp = df_reg.shape[0]
    block_size = np.ceil(n_snp / n_block).astype(int)

    # Places to make a split
    cut_set = set([0, n_snp])
    temp_v = df_reg["CHR"].values
    cut_set.update(np.arange(1, n_snp)[temp_v[1:] != temp_v[:-1]])

    # Places not to make a split
    nocut_set = set()
    for df_pannot in pannot_list:
        pAN = [x for x in df_pannot if x.startswith("pAN")][0]
        temp_dic = {x: y for x, y in zip(df_pannot["SNP"], df_pannot[pAN])}
        temp_v = np.array(
            [temp_dic[x] if x in temp_dic else sym_non_pAN for x in df_reg["SNP"]]
        )
        ind_select = (temp_v[1:] == temp_v[:-1]) & (temp_v[1:] != sym_non_pAN)
        nocut_set.update(np.arange(1, n_snp)[ind_select])

    dic_block, i_block, ind_s = {}, 0, 0
    for i in range(1, n_snp + 1):
        if (i in cut_set) | ((i - ind_s > block_size) & (i not in nocut_set)):
            dic_block[i_block] = (ind_s, i)
            ind_s = i
            i_block += 1
    return dic_block


def regress(
    df_reg,
    dic_block,
    n_sample_zsq,
    verbose=False,
    verbose_prefix='',
):

    """
    Single-pass regression (being called by the multi-stage regression `regress`).

    Parameters
    ----------
    df_reg : pd.DataFrame
        GDReg LD and DLD scores, with columns ['CHR', 'SNP', 'BP', 'ZSQ',
        'LD:AN1', 'LD:AN2', 'DLD:PAN:AN1', 'DLD:PAN:AN2'].
    dic_block : dict
        Block information. `dic_block[i] = (ind_s, ind_e)`.
    n_sample_zsq : int
        Number of samples used to compute the sumstats.
    dic_block : dict
        Block information. `dic_block[i] = [ind_s, ind_e]`.
    flag_tau_only : bool
        If true, only regress for \tau (not \tau and \rho).
    verbose : bool, default=False
        If to output messages.

    Returns
    -------
    dic_res_reg : dict
        Regression results.

        - dic_res_reg['term'] : list of terms.
        - dic_res_reg['coef'] : estimated coefs. np.ndarray(dtype=np.float32).
        - dic_res_reg['coef_jn'] : JN estimated coefs. np.ndarray(dtype=np.float32).
        - dic_res_reg['coef_jn_cov'] : estimated coef covariance. np.ndarray(dtype=np.float32).

    """

    start_time = time.time()

    n_snp = df_reg.shape[0]
    CHR_list = sorted(set(df_reg["CHR"]))

    LD_list = [x for x in df_reg if x.startswith("LD:")]
    DLD_list = [x for x in df_reg if x.startswith("DLD:")]
    reg_list = LD_list + DLD_list + ["E"]

    if verbose:
        print(verbose_prefix+"# Call: gdreg.regress.regress")
        print(
            verbose_prefix+
            "    n_snp=%d, n_block=%d, n_sample_zsq=%d"
            % (df_reg.shape[0], len(dic_block), n_sample_zsq)
        )
        print(verbose_prefix+"    %d regressors : %s" % (len(reg_list), ", ".join(reg_list)))

    # Regression weights
    # 1. LD : 1 / l_j
    # 2. Zsq variance : 1 / (1 + N h_g^2 l_j / M) ^ 2
    v_ld = df_reg["LD:AN:ALL"].values.clip(min=1)
    v_zsq_var = (
        n_sample_zsq * 0.5 * df_reg["LD:AN:ALL"].values / df_reg.shape[0]
        + 0.5 * df_reg["E"].values
    ).clip(min=0.1)
    v_w = np.sqrt(1 / v_ld / v_zsq_var)
    v_w = v_w * 0 + 1
    v_w = v_w / v_w.mean()

    # Regression
    mat_X = df_reg[reg_list].values.astype(np.float32)
    mat_X[:, :-1] *= n_sample_zsq
    v_y = df_reg["ZSQ"].values.astype(np.float32)

    mat_X = (mat_X.T * v_w).T
    v_y = v_y * v_w

    coef, coef_mean, coef_cov = reg_bjn(v_y, mat_X, dic_block)
    dic_res_reg = {
        "term": reg_list,
        "coef": coef,
        "coef_jn": coef_mean,
        "coef_jn_cov": coef_cov,
    }

    if verbose:
        print(verbose_prefix+"    Completed, time=%0.1fs" % (time.time() - start_time))
    return dic_res_reg


def reg_bjn(v_y, mat_X, dic_block, verbose=False):
    """
    OLS with block jackknife.

    Parameters
    ----------
    v_y : np.ndarray(dtype=np.float32)
        Response variable of shape (n_sample,).
    mat_X : np.ndarray(dtype=np.float32)
        Regressors of shape (n_sample, n_regressor).
    dic_block : dict
        Block information. `dic_block[i] = (ind_s, ind_e)`.
    verbose :  bool, default=False
        If to output messages.

    Returns
    -------
    coef_full : np.ndarray(dtype=np.float32)
        Estimates using full data of shape (n_regressor,).
    coef_mean : np.ndarray(dtype=np.float32)
        Jackknife bias-corrected estimates of shape (n_regressor,).
    coef_cov : np.ndarray(dtype=np.float32)
        Jackknife covariance of shape (n_regressor, n_regressor).
    """

    v_y = v_y.reshape([-1, 1])
    n_sample, n_regressor = mat_X.shape

    block_list = sorted(dic_block)
    n_block = len(block_list)
    v_block_size = [dic_block[x][1] - dic_block[x][0] for x in block_list]
    v_block_size = np.array(v_block_size, dtype=np.float32)

    # Check if dic_block covers all samples
    ind_select = np.zeros(n_sample, dtype=int)
    for block in block_list:
        ind_select[dic_block[block][0] : dic_block[block][1]] += 1

    n_0, n_g1 = (ind_select == 0).sum(), (ind_select > 1).sum()
    err_msg = "%d/%d of %d samples in no / >1 blocks" % (n_0, n_g1, n_sample)
    assert (n_0 == 0) & (n_g1 == 0), err_msg

    if verbose:
        print("# Call: gdreg.regress.reg_block_jn")
        print(
            "    n_sample=%d, n_regressor=%d, n_block=%d"
            % (n_sample, n_regressor, n_block)
        )

    # Full regression
    mat_xtx = np.dot(mat_X.T, mat_X) / n_sample
    mat_xty = np.dot(mat_X.T, v_y) / n_sample
    coef = np.linalg.solve(mat_xtx, mat_xty).reshape([-1])

    coef_block = np.zeros([n_block, n_regressor], dtype=np.float32)
    for i, block in enumerate(block_list):
        ind_s, ind_e = dic_block[block]
        mat_X_block = mat_X[ind_s:ind_e, :]
        v_y_block = v_y[ind_s:ind_e]

        mat_xtx_block = mat_xtx - np.dot(mat_X_block.T, mat_X_block) / n_sample
        mat_xty_block = mat_xty - np.dot(mat_X_block.T, v_y_block) / n_sample
        coef_block[i, :] = np.linalg.solve(mat_xtx_block, mat_xty_block).reshape([-1])
    #     print(coef_block)

    # Jacknife : mean & covariance
    v_h = n_sample / v_block_size
    coef_mean = (-coef_block + coef).sum(axis=0) + (coef_block.T / v_h).sum(axis=1)

    mat_tau = np.zeros([n_block, n_regressor], dtype=np.float32)
    for i in np.arange(n_block):
        mat_tau[i, :] = v_h[i] * coef - (v_h[i] - 1) * coef_block[i, :]

    coef_cov = np.zeros([n_regressor, n_regressor], dtype=np.float32)
    for i in np.arange(n_block):
        temp_v = mat_tau[i, :] - coef_mean
        coef_cov += np.outer(temp_v, temp_v) / n_block / (v_h[i] - 1)

    return coef, coef_mean, coef_cov
