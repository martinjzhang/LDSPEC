import numpy as np
import scipy as sp
import pandas as pd
import time
import gdreg
import warnings


def estimate(
    dic_data,
    df_score,
    df_sumstats,
    dic_annot_path={},
    dic_pannot_path={},
    flag_cross_term=False,
    n_jn_block=100,
    verbose=False,
):
    """
    GDREG estimation

    Parameters
    ----------
    dic_data : dict
        Genotype data reader, organized, for each CHR, as

        - dic_data[CHR]['pgen'] : .pgen file path
        - dic_data[CHR]['pvar'] : .pvar pd.DataFrame
        - dic_data[CHR]['psam'] : .psam pd.DataFrame
    df_score : pd.DataFrame, default=None
        GDREG LD and DLD scores, with columns ['CHR', 'SNP', 'BP', 'LD:AN:name1', 'LD:AN:name2',
        'LD:E', 'DLD:pAN:name1', 'DLD:pAN:name2', ...]. Must contain 'LD:AN:allXX' where
        XX is one of ["", "_common", "_ld", "_rare"].
    df_sumstats : pd.DataFrame
        Summary statistics with columns ['SNP', 'N', 'Z', 'A1', 'A2']
    dic_annot_path : dic of dic of strs
        File path for single-SNP annotation. dic_annot_path[annot_name][CHR] contains the
        `.annot.gz` file path for annot_file `annot_name` and CHR `CHR`.
    dic_pannot_path : dic of dic of strs
        File path for SNP-pair annotation. dic_pannot_path[annot_name][CHR] contains the
        `.pannot_mat.npz` file path for annotation pAN and and CHR `CHR`. Dimension of the
        sparse matrix should match `dic_data[CHR][pvar]`.
    flag_cross_term : bool, default=False
        If True, also use cross terms (Z_i Z_j) for regression, for SNP pairs i,j within
        10000 SNPs and covered by at least one pannot.
    n_jn_block : int, default=100
        Number of JN blocks.

    Returns
    -------
    df_res : pd.DataFrame
        GDREG results.

        - dic_res[0] : regression using only LD scores.
            - dic_res['term'] : regression terms
            - dic_res['coef'] : regression coefficients
            - dic_res['coef_jn'] : JN-debiased regression coefficients
            - dic_res['coef_jn_cov'] : regression coefficient covariance based on JN
            - dic_res['summary'] : result summary
        - dic_res[1] : regression using both LD and DLD scores, structure same as dic_res[0]


    TODO
    ----
    - Only do the joint analyses
    """

    start_time = time.time()
    if verbose:
        print("# Call: gdreg.regress.estimate")

    # Annotations
    CHR_list = sorted(dic_data)  # CHR_list contains all CHRs
    for annot_name in dic_annot_path:
        err_msg = "%s does not contain all CHRs in dic_data" % annot_name
        assert set(dic_annot_path[annot_name]) == set(CHR_list), err_msg
    for annot_name in dic_pannot_path:
        err_msg = "%s does not contain all CHRs in dic_data" % annot_name
        assert set(dic_pannot_path[annot_name]) == set(CHR_list), err_msg

    # df_score
    df_score.index = df_score["SNP"]
    ind_rm = df_score.isna().sum(axis=1) > 0
    df_score = df_score.loc[~ind_rm].copy()
    LD_list = [x for x in df_score if x.startswith("LD:")]
    DLD_list = [x for x in df_score if x.startswith("DLD:")]
    assert "E" in df_score, "'E' not in df_score"
    if verbose:
        print(
            "    df_score : remove %d rows with NA values, %d remaining"
            % (ind_rm.sum(), df_score.shape[0])
        )
        print("        %d LD scores, %s DLD scores" % (len(LD_list), len(DLD_list)))

    # df_sumstats
    n_sample_zsq = df_sumstats["N"].mean().astype(int)
    dic_zsc = {x: y for x, y in zip(df_sumstats["SNP"], df_sumstats["Z"])}
    outlier_thres = max(80, 0.001 * n_sample_zsq)  # Finucane 2015 Nat Genet
    dic_zsc = {x: y for x, y in dic_zsc.items() if y ** 2 < outlier_thres}

    if verbose:
        print(
            "    df_sumstats : n_snp=%d, n_sample_zsq=%d"
            % (df_sumstats.shape[0], n_sample_zsq)
        )
        print(
            "        Remove duplicate or ZSQ>%0.1f SNPs, %d remaining, avg. zsq=%0.2f"
            % (
                outlier_thres,
                len(dic_zsc),
                np.mean(np.array(list(dic_zsc.values())) ** 2),
            )
        )

    # df_reg
    if flag_cross_term:
        df_reg = df_score[["SNP", "CHR", "BP"]].copy()
    else:
        ind_select = ["|" not in x for x in df_score["SNP"]]
        df_reg = df_score.loc[ind_select, ["SNP", "CHR", "BP"]].copy()
    df_reg.drop_duplicates("SNP", inplace=True)
    df_reg["SNP1"] = [x.split("|")[0] for x in df_reg["SNP"]]
    df_reg["SNP2"] = [x.split("|")[-1] for x in df_reg["SNP"]]
    df_reg = df_reg.loc[(df_reg["SNP1"].isin(dic_zsc)) & (df_reg["SNP2"].isin(dic_zsc))]
    df_reg["ZSQ"] = [
        dic_zsc[x] * dic_zsc[y] for x, y in zip(df_reg["SNP1"], df_reg["SNP2"])
    ]
    df_reg.index = df_reg["SNP"]
    df_reg.sort_values(by=["CHR", "BP"], inplace=True)

    if verbose:
        print("    Regression : n_rows=%d, n_block=%d" % (df_reg.shape[0], n_jn_block))

    # Regression
    dic_res = {}
    temp_df_reg = df_reg.join(df_score[LD_list + DLD_list + ["E"]])
    dic_block = get_block(temp_df_reg, n_block=n_jn_block)
    dic_res = regress(
        temp_df_reg,
        dic_block,
        n_sample_zsq,
        verbose=verbose,
        verbose_prefix="    ",
    )
    dic_res["summary"] = summarize(
        dic_res,
        dic_data,
        dic_annot_path=dic_annot_path,
        dic_pannot_path=dic_pannot_path,
    )

    if verbose:
        print("    Completed, time=%0.1fs" % (time.time() - start_time))

    return dic_res


def summarize(
    dic_res,
    dic_data,
    dic_annot_path={},
    dic_pannot_path={},
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
    dic_data : dict
        Genotype data reader. Must contain SNPs from all chromosomes.

        - dic_data[CHR]['pgen'] : .pgen file path
        - dic_data[CHR]['pvar'] : .pvar pd.DataFrame
        - dic_data[CHR]['psam'] : .psam pd.DataFrame
    dic_annot_path : dic of dic of strs
        File path for single-SNP annotation. dic_annot_path[annot_name][CHR] contains the
        `.annot.gz` file path for annot_file `annot_name` and CHR `CHR`.
    dic_pannot_path : dic of dic of strs
        File path for SNP-pair annotation. dic_pannot_path[annot_name][CHR] contains the
        `.pannot_mat.npz` file path for annotation pAN and and CHR `CHR`. Dimension of the
        sparse matrix should match `dic_data[CHR][pvar]`.

    Returns
    -------
    df_summary : pd.DataFrame
        Regression result summary.

        - tau,tau_se : \tau coefficient
        - h2,h2_se : total heritability of SNPs in a given annot
        - enrich,enrich_se : heritability enrichment
        - rho,rho_se : \rho coefficient
        - cov,cov_se : total covariance of SNP pairs in a given pannot
        - r2,r2_se : `cov` divided by total \sqrt{ h_ps_i h_ps_j } of SNP pairs in a given pannot,
            where h_ps_i is the per-SNP heritablity of SNP i. `r2_se` based on if `cov` is
            significantly different from 0.

    TODO
    ----

    """

    # Annots and pannots
    CHR_list = sorted(dic_data)  # CHR_list contains all CHRs
    for annot_name in dic_annot_path:
        err_msg = "%s does not contain all CHRs in dic_data" % annot_name
        assert set(dic_annot_path[annot_name]) == set(CHR_list), err_msg
    for annot_name in dic_pannot_path:
        err_msg = "%s does not contain all CHRs in dic_data" % annot_name
        assert set(dic_pannot_path[annot_name]) == set(CHR_list), err_msg
        
    AN_list = []
    for annot_name in dic_annot_path:
        CHR = list(dic_annot_path[annot_name])[0]
        temp_df = gdreg.util.read_annot(dic_annot_path[annot_name][CHR], nrows=5)
        AN_list.extend([x for x in temp_df if x.startswith("AN:")])
    pAN_list = list(dic_pannot_path)

    # Results info
    res_LD_list = [x for x in dic_res["term"] if x.startswith("LD:")]
    res_AN_list = [x.replace("LD:", "") for x in res_LD_list]
    res_DLD_list = [x for x in dic_res["term"] if x.startswith("DLD:")]
    res_pAN_list = sorted(
        set([x.replace("DLD:", "").split("|")[0] for x in res_DLD_list])
    )

    # Check consistency between results and annotation df_annot
    err_msg = "df_annot does not contain all annots in dic_res"
    assert len(set(res_AN_list) - set(AN_list)) == 0, err_msg
    err_msg = "dic_pannot_mat does not contain all pannots in dic_res"
    assert len(set(res_pAN_list) - set(pAN_list)) == 0, err_msg

    # Summary
    dic_coef = {x: y for x, y in zip(dic_res["term"], dic_res["coef_jn"])}
    v_se = np.sqrt(np.diag(dic_res["coef_jn_cov"]))
    dic_coef_se = {x: y for x, y in zip(dic_res["term"], v_se)}
    df_sum_tau = pd.DataFrame(
        index=res_AN_list,
        data={
            "annot": res_AN_list,
            "n_snp": 0,
            "tau": [dic_coef["LD:%s" % x] for x in res_AN_list],
            "tau_se": [dic_coef_se["LD:%s" % x] for x in res_AN_list],
            "h2": np.nan,
            "h2_se": np.nan,
            "enrich": np.nan,
            "enrich_se": np.nan,
        },
    )
    df_sum_rho = pd.DataFrame(
        index=res_pAN_list,
        data={
            "pannot": res_pAN_list,
            "n_pair": 0,
            "rho": [dic_coef["DLD:%s" % x] for x in res_pAN_list],
            "rho_se": [dic_coef_se["DLD:%s" % x] for x in res_pAN_list],
            "cov": np.nan,  # Avg. cov
            "cov_se": np.nan,
            "r2": np.nan,  # Avg. cor
            "r2_se": np.nan,
        },
    )
    
#     return {"tau": df_sum_tau, "rho": df_sum_rho}

    # Iterate over CHR_list to collect info
    dic_AN_n_snp = {x: 0 for x in res_AN_list}
    dic_AN_n_snp_ref = {x: 0 for x in res_AN_list}
    dic_AN_type = {x: "binary" for x in res_AN_list}
    dic_AN_v = {x: np.zeros(len(res_AN_list), dtype=np.float32) for x in res_AN_list}
    dic_AN_v_ref = {
        x: np.zeros(len(res_AN_list), dtype=np.float32) for x in res_AN_list
    }

    dic_pAN_n_pair = {x: 0 for x in res_pAN_list}
    dic_pAN_v = {x: np.zeros(len(res_pAN_list), dtype=np.float32) for x in res_pAN_list}
    dic_pAN_var_total = {x: 0 for x in res_pAN_list}

    for CHR in CHR_list:
        # df_annot_chr
        df_annot_chr = dic_data[CHR]["pvar"][["SNP"]].copy()
        for annot_name in dic_annot_path:
            temp_df = gdreg.util.read_annot(dic_annot_path[annot_name][CHR])
            for AN in [x for x in temp_df if x.startswith("AN:")]:
                if AN in res_AN_list:
                    temp_dic = {
                        x: y for x, y in zip(temp_df["SNP"], temp_df[AN]) if y != 0
                    }
                    df_annot_chr[AN] = [
                        temp_dic[x] if x in temp_dic else 0 for x in df_annot_chr["SNP"]
                    ]
                    df_annot_chr[AN] = df_annot_chr[AN].astype(np.float32)
        # Update info
        for AN in res_AN_list:
            dic_AN_n_snp[AN] = dic_AN_n_snp[AN] + (df_annot_chr[AN] == 1).sum()
            if len(df_annot_chr[AN].unique()) > 2:
                dic_AN_type[AN] = "non-binary"
            dic_AN_v[AN] = (
                dic_AN_v[AN]
                + df_annot_chr.loc[df_annot_chr[AN] == 1, res_AN_list]
                .sum(axis=0)
                .values
            )

            ref_col_list = res_AN_list
            for term in ["_common", "_lf", "_rare"]:  # if not, use all SNPs as ref
                if AN.endswith(term):
                    ref_col_list = [x for x in res_AN_list if x.endswith(term)]
            ind_ref = (df_annot_chr[ref_col_list].values == 1).sum(axis=1) > 0
            dic_AN_n_snp_ref[AN] = dic_AN_n_snp_ref[AN] + ind_ref.sum()
            dic_AN_v_ref[AN] = (
                dic_AN_v_ref[AN]
                + df_annot_chr.loc[ind_ref, res_AN_list].sum(axis=0).values
            )

        # dic_mat_G_chr
        dic_mat_G_chr = {}
        for pAN in res_pAN_list:
            dic_mat_G_chr[pAN] = gdreg.util.read_pannot_mat(dic_pannot_path[pAN][CHR])

        # Update info
        v_persnp_h2 = np.zeros(dic_data[CHR]["pvar"].shape[0], dtype=np.float32)
        for AN in res_AN_list:
            v_persnp_h2 = v_persnp_h2 + df_sum_tau.loc[AN, "tau"] * df_annot_chr[AN]
        v_persnp_h2_sqrt = np.sqrt(v_persnp_h2).astype(np.float32)

        for pAN in res_pAN_list:
            dic_pAN_n_pair[pAN] = dic_pAN_n_pair[pAN] + dic_mat_G_chr[pAN].sum()
            temp_list = [
                (dic_mat_G_chr[pAN] * dic_mat_G_chr[x]).sum() for x in res_pAN_list
            ]
            dic_pAN_v[pAN] = dic_pAN_v[pAN] + np.array(temp_list, dtype=np.float32)
            dic_pAN_var_total[pAN] = dic_pAN_var_total[pAN] + dic_mat_G_chr[pAN].dot(
                v_persnp_h2_sqrt
            ).T.dot(v_persnp_h2_sqrt)

    # Summary : n_snp,n_pair
    df_sum_tau["n_snp"] = [dic_AN_n_snp[x] for x in res_AN_list]
    df_sum_rho["n_pair"] = [dic_pAN_n_pair[x] for x in res_pAN_list]

    # Summary : h2, h2_se, enrich, enrich_se
    df_cov = pd.DataFrame(
        index=dic_res["term"], columns=dic_res["term"], data=dic_res["coef_jn_cov"]
    )
    temp_mat = df_cov.loc[res_LD_list, res_LD_list].values
    for AN in res_AN_list:
        if dic_AN_type[AN] != "binary":
            continue

        n_snp_AN = dic_AN_n_snp[AN]
        n_snp_ref = dic_AN_n_snp_ref[AN]

        # h2, h2_se
        df_sum_tau.loc[AN, "h2"] = (dic_AN_v[AN] * df_sum_tau["tau"]).sum()
        df_sum_tau.loc[AN, "h2_se"] = np.sqrt(dic_AN_v[AN].dot(temp_mat).dot(dic_AN_v[AN]))

        # enrich, enrich_se
        n_snp_dif = n_snp_ref - n_snp_AN
        if n_snp_dif < n_snp_AN * 0.1:
            continue

        temp_v_combine = (
            dic_AN_v[AN] * (1 / n_snp_AN + 1 / n_snp_dif) - dic_AN_v_ref[AN] / n_snp_dif
        )

        h2_ps_ref = (dic_AN_v_ref[AN] * df_sum_tau["tau"]).sum() / n_snp_ref
        df_sum_tau.loc[AN, "enrich"] = df_sum_tau.loc[AN, "h2"] / n_snp_AN / h2_ps_ref

        dif_ = (temp_v_combine * df_sum_tau["tau"]).sum()
        se_ = np.sqrt(temp_v_combine.dot(temp_mat).dot(temp_v_combine))
        df_sum_tau.loc[AN, "enrich_se"] = np.absolute(
            se_ / dif_ * (df_sum_tau.loc[AN, "enrich"] - 1)
        )

    # Summary : cov, cov_se, r2, r2_se
    temp_mat = df_cov.loc[res_DLD_list, res_DLD_list].values
    for pAN in res_pAN_list:
        # cov, cov_se
        df_sum_rho.loc[pAN, "cov"] = (dic_pAN_v[pAN] * df_sum_rho["rho"]).sum()
        df_sum_rho.loc[pAN, "cov_se"] = np.sqrt(dic_pAN_v[pAN].dot(temp_mat).dot(dic_pAN_v[pAN]))

        # r2 and r2_se
        df_sum_rho.loc[pAN, "r2"] = df_sum_rho.loc[pAN, "cov"] / dic_pAN_var_total[pAN]
        df_sum_rho.loc[pAN, "r2_se"] = df_sum_rho.loc[pAN, "cov_se"] / dic_pAN_var_total[pAN]

    return {"tau": df_sum_tau, "rho": df_sum_rho}


def get_block(df_reg, pannot_list=[], sym_non_pAN="non-pAN", n_block=100):
    """
    Block partition. SNPs on the same .pannot annotation go to the same block.
    SNPs on different CHRs go to different blocks.

    Parameters
    ----------
    df_reg : pd.DataFrame
        SNPs used in regression. Must be sorted by genomic location. Must contain
        ['CHR', 'SNP', 'SNP1', 'SNP2'].
    pannot_list : list of pd.DataFrame, default=[]
        Each element corresponds to a SNP-pair annotation with columns
        ['CHR', 'SNP', 'BP', 'pAN:name']
    sym_non_pAN : str, default='non-pAN'
        Symbol for SNPs not in the SNP-pair annotation.
    n_block: int, default=100
        Number of jackknife blocks.

    Returns
    -------
    dic_block : dict
        Block information. `dic_block[i] = (ind_s, ind_e)`.

    TODO
    ----
    1. Currently based only on pannot and 'SNP1'.
    """

    n_snp = df_reg.shape[0]
    block_size = np.ceil(n_snp / n_block).astype(int)

    # Split locations
    cut_set = set([0, n_snp])
    temp_v = df_reg["CHR"].values
    cut_set.update(np.arange(1, n_snp)[temp_v[1:] != temp_v[:-1]])

    # Non-split locations
    nocut_set = set()
    for df_pannot in pannot_list:
        pAN = [x for x in df_pannot if x.startswith("pAN")][0]
        temp_dic = {x: y for x, y in zip(df_pannot["SNP"], df_pannot[pAN])}

        # SNP1
        temp_v = np.array(
            [temp_dic[x] if x in temp_dic else sym_non_pAN for x in df_reg["SNP1"]]
        )
        ind_select1 = (temp_v[1:] == temp_v[:-1]) & (temp_v[1:] != sym_non_pAN)
        # SNP2
        temp_v = np.array(
            [temp_dic[x] if x in temp_dic else sym_non_pAN for x in df_reg["SNP2"]]
        )
        ind_select2 = (temp_v[1:] == temp_v[:-1]) & (temp_v[1:] != sym_non_pAN)

        nocut_set.update(np.arange(1, n_snp)[ind_select1 | ind_select2])

    #         temp_v = np.array(
    #             [temp_dic[x] if x in temp_dic else sym_non_pAN for x in df_reg["SNP1"]]
    #         )
    #         ind_select = (temp_v[1:] == temp_v[:-1]) & (temp_v[1:] != sym_non_pAN)
    #         nocut_set.update(np.arange(1, n_snp)[ind_select])

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
    verbose_prefix="",
):

    """
    GDREG regression

    Parameters
    ----------
    df_reg : pd.DataFrame
        GDReg LD and DLD scores, with columns ['CHR', 'SNP', 'BP', 'SNP1', 'SNP2', 'ZSQ',
        'LD:AN:name1', 'LD:AN:name2', 'DLD:pAN:name1', 'DLD:pAN:name2', ...].
        If `SNP1|SNP2` in `df_reg`, `SNP1` and `SNP2` must both be in `df_reg`.
    dic_block : dict
        Block information. `dic_block[i] = (ind_s, ind_e)`.
    n_sample_zsq : int
        Number of samples used to compute the sumstats.

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

    #     n_snp = df_reg.shape[0]
    #     CHR_list = sorted(set(df_reg["CHR"]))

    LD_list = [x for x in df_reg if x.startswith("LD:")]
    DLD_list = [x for x in df_reg if x.startswith("DLD:")]
    reg_list = LD_list + DLD_list + ["E"]

    if verbose:
        print(verbose_prefix + "# Call: gdreg.regress.regress")
        print(
            verbose_prefix
            + "    n_snp=%d, n_block=%d, n_sample_zsq=%d"
            % (df_reg.shape[0], len(dic_block), n_sample_zsq)
        )
        print(
            verbose_prefix
            + "    %d regressors : %s" % (len(reg_list), ", ".join(reg_list))
        )

    # Regression weights
    # TODO : Balance Z_i^2 and Z_i Z_j
    # 1. LD weights (clipped at max=10)
    #    - Z_i^2 : 1 / l_i, where l_i = \sum_j r_ij^2
    #    - Z_i Z_j : 1 / l_ij, where l_ij = \sum_k r_ik r_jk
    temp_list = [
        x for x in df_reg if x.startswith("LD:AN:all") | x.startswith("LD:AN:ALL")
    ]
    v_ld = df_reg[temp_list].sum(axis=1).values.clip(min=0.1)
    # 2. Zsq variance weights : (clipped at max=10)
    #    - Z_i^2 : 1 / [ 2 (N l_i / M + 1) ^ 2 ]
    #    - Z_i Z_j : 1 / [ (N l_i / M + 1) (N l_j / M + 1) + (N l_ij / M + 1) ^ 2 ]
    #
    n_snp = (df_reg["SNP1"] == df_reg["SNP2"]).sum()
    dic_ld = {(x, y): z for x, y, z in zip(df_reg["SNP1"], df_reg["SNP2"], v_ld)}
    v_zsq_var = [
        (n_sample_zsq * dic_ld[(s1, s1)] / n_snp + 1)
        * (n_sample_zsq * dic_ld[(s2, s2)] / n_snp + 1)
        + (n_sample_zsq * dic_ld[(s1, s2)] / n_snp + ld_e) ** 2
        for s1, s2, ld_e in zip(df_reg["SNP1"], df_reg["SNP2"], df_reg["E"])
    ]
    v_zsq_var = np.array(v_zsq_var, dtype=np.float32).clip(min=0.1)

    v_w = np.sqrt(1 / v_ld / v_zsq_var).astype(np.float32)
    v_w = v_w / v_w.mean()

    #     # Regression weights
    #     # 1. LD : 1 / l_j
    #     # 2. Zsq variance : 1 / (1 + N h_g^2 l_j / M) ^ 2
    #     LD_all_list = [
    #         x for x in df_reg if x.startswith("LD:AN:all") | x.startswith("LD:AN:ALL")
    #     ]
    #     v_ld = df_reg[LD_all_list].sum(axis=1).values.clip(min=1)
    #     v_zsq_var = (
    #         n_sample_zsq * 0.5 * df_reg[LD_all_list].sum(axis=1).values / df_reg.shape[0]
    #         + 0.5 * df_reg["E"].values
    #     ).clip(min=0.1)
    #     v_w = np.sqrt(1 / v_ld / v_zsq_var)
    #     v_w = v_w / v_w.mean()

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
        print(
            verbose_prefix + "    Completed, time=%0.1fs" % (time.time() - start_time)
        )
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
