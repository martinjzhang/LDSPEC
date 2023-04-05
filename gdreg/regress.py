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
    dic_avgr={},
    flag_cross_term=False,
    flag_nofil_snp=False,
    n_jn_block=100,
    null_model=None,
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
        XX is one of ["", "_common", "_ld", "_rare"]. pannots with 'prox' are used as reference
        for computing ecov and ecor
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
    flag_nofil_snp : bool, default=False
        If True, turning off outlier filter.
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
    - Add jn results for other quantities
    """

    start_time = time.time()
    if verbose:
        print("# Call: gdreg.regress.estimate")

    # Annotations
    CHR_list = sorted(dic_data)  # CHR_list contains all CHRs
    n_jn_block = min(n_jn_block, len(CHR_list) * 10)  # <=10 JN blocks for each CHR
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
    if flag_nofil_snp:
        outlier_thres = 1e8
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
        null_model=null_model,
        verbose=verbose,
        verbose_prefix="    ",
    )

    # Summary only in the estimation mode
    if null_model is None:
        dic_res["summary"] = summarize(
            dic_res,
            dic_data,
            dic_annot_path=dic_annot_path,
            dic_pannot_path=dic_pannot_path,
            dic_avgr=dic_avgr,
            verbose=verbose,
            verbose_prefix="    ",
        )

    if verbose:
        print("    Completed, time=%0.1fs" % (time.time() - start_time))

    return dic_res


def summarize(
    dic_res,
    dic_data,
    dic_annot_path={},
    dic_pannot_path={},
    dic_avgr={},
    verbose=False,
    verbose_prefix="",
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
    dic_annot_path : dic of dic of strs, , default={}
        File path for single-SNP annotation. dic_annot_path[annot_name][CHR] contains the
        `.annot.gz` file path for annot_file `annot_name` and CHR `CHR`.
    dic_pannot_path : dic of dic of strs, default={}
        File path for SNP-pair annotation. dic_pannot_path[annot_name][CHR] contains the
        `.pannot_mat.npz` file path for annotation pAN and and CHR `CHR`. Dimension of the
        sparse matrix should match `dic_data[CHR][pvar]`.
    dic_avgr : dic, default={}
        dic_avgr[pAN] contains the average LD across all pairs in pAN.

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

    start_time = time.time()

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

    if len(dic_avgr) == 0:
        dic_avgr = {x: 0 for x in pAN_list}

    # dis_res
    res_AN_list = [x.replace("LD:", "") for x in dic_res["term"] if x.startswith("LD:")]
    res_pAN_list = [
        x.replace("DLD:", "") for x in dic_res["term"] if x.startswith("DLD:")
    ]
    res_prox_list = [x for x in res_pAN_list if "prox" in x]  # pAN's for priximity
    err_msg = "df_annot does not contain all annots in dic_res"
    assert len(set(res_AN_list) - set(AN_list)) == 0, err_msg
    err_msg = "dic_pannot_mat does not contain all pannots in dic_res"
    assert len(set(res_pAN_list) - set(pAN_list)) == 0, err_msg
    # TODO : check dic_avgr

    if verbose:
        print(verbose_prefix + "# Call: gdreg.regress.summarize")
        print(
            verbose_prefix
            + "    %d annots, %d pannots" % (len(res_AN_list), len(res_pAN_list))
        )
        print(
            verbose_prefix
            + "    %d prox pannots : %s"
            % (len(res_prox_list), ", ".join(res_prox_list))
        )

    temp_list = [x.replace("DLD:", "").replace("LD:", "") for x in dic_res["term"]]
    dic_coef = {x: y for x, y in zip(temp_list, dic_res["coef"])}
    dic_coef_jn = {x: y for x, y in zip(temp_list, dic_res["coef_jn"])}
    df_coef_cov = pd.DataFrame(
        index=temp_list,
        columns=temp_list,
        data=dic_res["coef_jn_cov"],
        dtype=np.float32,
    )
    n_jn_block = dic_res["coef_block"].shape[0]
    df_coef_block = pd.DataFrame(
        columns=temp_list, data=dic_res["coef_block"], dtype=np.float32
    )  # (n_jn_block, n_coef)

    # Summary : dic_res and dfs
    df_sum_tau = pd.DataFrame(
        index=res_AN_list,
        data={
            "annot": res_AN_list,
            "type": "",
            "n_snp": 0,
            "tau": [dic_coef_jn[x] for x in res_AN_list],
            "tau_se": [np.sqrt(df_coef_cov.loc[x, x]) for x in res_AN_list],
            "tau_p": np.nan,
            "h2": np.nan,  # Total heritability
            "h2_se": np.nan,
            "h2_p": np.nan,
            "h2s": np.nan,  # Sum of causal effect size variance (SCV)
            "h2s_se": np.nan,
            "h2s_p": np.nan,
            "h2_enrich": np.nan,
            "h2_enrich_se": np.nan,
            "h2_enrich_p": np.nan,
            "h2s_enrich": np.nan,
            "h2s_enrich_se": np.nan,
            "h2s_enrich_p": np.nan,
            "h2_shrink": np.nan,
            "h2_shrink_se": np.nan,
            "h2_shrink_p": np.nan,
        },
    )

    df_sum_rho = pd.DataFrame(
        index=res_pAN_list,
        data={
            "pannot": res_pAN_list,
            "n_pair": 0,
            "rho": [dic_coef_jn[x] for x in res_pAN_list],
            "rho_se": [np.sqrt(df_coef_cov.loc[x, x]) for x in res_pAN_list],
            "rho_p": np.nan,
            "cov": np.nan,  # Total covariance
            "cov_se": np.nan,
            "cov_p": np.nan,
            "cor": np.nan,  # Average correlation
            "cor_se": np.nan,
            "cor_p": np.nan,
            "ecov": np.nan,  # Excess total covariance
            "ecov_se": np.nan,
            "ecov_p": np.nan,
            "ecor": np.nan,  # Excess average correlation
            "ecor_se": np.nan,
            "ecor_p": np.nan,
        },
    )

    dic_jn = {
        "v_h": dic_res["v_h"],
        "res_AN_list": res_AN_list,
        "h2": np.zeros(len(res_AN_list), dtype=np.float32),
        "h2.jn": np.zeros([n_jn_block, len(res_AN_list)], dtype=np.float32),
        "h2_enrich": np.zeros(len(res_AN_list), dtype=np.float32),
        "h2_enrich.jn": np.zeros([n_jn_block, len(res_AN_list)], dtype=np.float32),
        "h2s": np.zeros(len(res_AN_list), dtype=np.float32),
        "h2s.jn": np.zeros([n_jn_block, len(res_AN_list)], dtype=np.float32),
        "h2s_enrich": np.zeros(len(res_AN_list), dtype=np.float32),
        "h2s_enrich.jn": np.zeros([n_jn_block, len(res_AN_list)], dtype=np.float32),
        "res_pAN_list": res_pAN_list,
        "cov": np.zeros(len(res_pAN_list), dtype=np.float32),
        "cov.jn": np.zeros([n_jn_block, len(res_pAN_list)], dtype=np.float32),
        "cor": np.zeros(len(res_pAN_list), dtype=np.float32),
        "cor.jn": np.zeros([n_jn_block, len(res_pAN_list)], dtype=np.float32),
        "ecov": np.zeros(len(res_pAN_list), dtype=np.float32),
        "ecov.jn": np.zeros([n_jn_block, len(res_pAN_list)], dtype=np.float32),
        "ecor": np.zeros(len(res_pAN_list), dtype=np.float32),
        "ecor.jn": np.zeros([n_jn_block, len(res_pAN_list)], dtype=np.float32),
    }

    # Iterate over CHR_list to collect info
    dic_AN_n_snp = {x: 0 for x in res_AN_list}
    dic_AN_type = {x: "binary" for x in res_AN_list}
    dic_AN_v = {
        x: np.zeros(len(res_AN_list), dtype=np.float32) for x in res_AN_list
    }  # n_overlap between all AN and the given AN
    dic_AN_n_snp_ref = {x: 0 for x in res_AN_list}
    dic_AN_v_ref = {
        x: np.zeros(len(res_AN_list), dtype=np.float32) for x in res_AN_list
    }  # All SNPs from the same common/lf/rare mbin as `AN`

    dic_AN_v_p = {
        x: np.zeros(len(res_pAN_list), dtype=np.float32) for x in res_AN_list
    }  # h2p coefficients
    dic_AN_v_p_ref = {
        x: np.zeros(len(res_pAN_list), dtype=np.float32) for x in res_AN_list
    }  # Reference h2p coefficients

    dic_pAN_n_pair = {x: 0 for x in res_pAN_list}
    dic_pAN_v = {
        x: np.zeros(len(res_pAN_list), dtype=np.float32) for x in res_pAN_list
    }  # n_overlap between all pAN and the given pAN
    dic_pAN_var = {x: 0 for x in res_pAN_list}  # Total sqrt(var_i var_j)
    dic_pAN_var_block = {x: [0] * n_jn_block for x in res_pAN_list}

    for CHR in CHR_list:
        # df_annot_chr
        df_annot_chr = pd.DataFrame(
            index=dic_data[CHR]["pvar"]["SNP"],
            columns=res_AN_list,
            data=0,
            dtype=np.float32,
        )
        for annot_name in dic_annot_path:
            temp_df = gdreg.util.read_annot(dic_annot_path[annot_name][CHR])
            for AN in [x for x in temp_df if x in res_AN_list]:
                temp_dic = {x: y for x, y in zip(temp_df["SNP"], temp_df[AN]) if y != 0}
                df_annot_chr[AN] = np.array(
                    [temp_dic[x] if x in temp_dic else 0 for x in df_annot_chr.index],
                    dtype=np.float32,
                )

        # Square-root per-SNP heritability : v_h2ps, v_h2ps_jn, mat_h2ps_block
        v_h2ps = np.zeros(df_annot_chr.shape[0], dtype=np.float32)
        v_h2ps_jn = np.zeros(df_annot_chr.shape[0], dtype=np.float32)
        mat_h2ps_block = np.zeros(
            [n_jn_block, df_annot_chr.shape[0]], dtype=np.float32
        )  # (n_jn_block, n_snp_chr)
        for AN in res_AN_list:
            v_h2ps += dic_coef[AN] * df_annot_chr[AN].values
            v_h2ps_jn += dic_coef_jn[AN] * df_annot_chr[AN].values
            for i in range(n_jn_block):
                mat_h2ps_block[i, :] += (
                    df_coef_block.loc[i, AN] * df_annot_chr[AN].values
                )
        v_h2ps = np.sqrt(v_h2ps.clip(min=0)).astype(np.float32)
        v_h2ps_jn = np.sqrt(v_h2ps_jn.clip(min=0)).astype(np.float32)
        mat_h2ps_block = np.sqrt(mat_h2ps_block.clip(min=0)).astype(np.float32)

        # Update annot info
        for AN in res_AN_list:
            dic_AN_n_snp[AN] += df_annot_chr[AN].sum()
            if len(df_annot_chr[AN].unique()) > 2:
                dic_AN_type[AN] = "non-binary"
            dic_AN_v[AN] += (
                df_annot_chr.loc[df_annot_chr[AN] == 1, res_AN_list].sum(axis=0).values
            )

            ref_col_list = res_AN_list
            for term in ["_common", "_lf", "_rare"]:  # if not, use all SNPs as ref
                if AN.endswith(term):
                    ref_col_list = [x for x in res_AN_list if x.endswith(term)]
            ind_ref = (df_annot_chr[ref_col_list].values == 1).sum(axis=1) > 0
            dic_AN_n_snp_ref[AN] += ind_ref.sum()
            dic_AN_v_ref[AN] += (
                df_annot_chr.loc[ind_ref, res_AN_list].sum(axis=0).values
            )

        # dic_mat_G_chr
        dic_mat_G_chr = {}
        for pAN in res_pAN_list:
            dic_mat_G_chr[pAN] = gdreg.util.read_pannot_mat(dic_pannot_path[pAN][CHR])

        # Update pannot info
        for pAN in res_pAN_list:
            dic_pAN_n_pair[pAN] += dic_mat_G_chr[pAN].sum()
            temp_list = [
                dic_mat_G_chr[pAN].multiply(dic_mat_G_chr[x]).sum()
                for x in res_pAN_list
            ]
            dic_pAN_v[pAN] += np.array(temp_list, dtype=np.float32)

            dic_pAN_var[pAN] += dic_mat_G_chr[pAN].dot(v_h2ps).T.dot(v_h2ps)
            for i in range(n_jn_block):
                dic_pAN_var_block[pAN][i] += (
                    dic_mat_G_chr[pAN]
                    .dot(mat_h2ps_block[i, :])
                    .T.dot(mat_h2ps_block[i, :])
                )

        # Update annot * pannot info for h2p and h2p_enrich
        for AN in res_AN_list:
            ref_col_list = res_AN_list
            for term in ["_common", "_lf", "_rare"]:  # if not, use all SNPs as ref
                if AN.endswith(term):
                    ref_col_list = [x for x in res_AN_list if x.endswith(term)]
            ind_ref = ((df_annot_chr[ref_col_list].values == 1).sum(axis=1) > 0) * 1
            for i_pAN, pAN in enumerate(res_pAN_list):
                dic_AN_v_p[AN][i_pAN] += dic_avgr[pAN] * (
                    dic_mat_G_chr[pAN]
                    .dot(df_annot_chr[AN].values)
                    .T.dot(df_annot_chr[AN].values)
                )
                dic_AN_v_p_ref[AN][i_pAN] += dic_avgr[pAN] * (
                    dic_mat_G_chr[pAN].dot(ind_ref).T.dot(ind_ref)
                )

    # Summary : n_snp,type,n_pair
    df_sum_tau["n_snp"] = [dic_AN_n_snp[x] for x in res_AN_list]
    df_sum_tau["type"] = [dic_AN_type[x] for x in res_AN_list]
    df_sum_rho["n_pair"] = [dic_pAN_n_pair[x] for x in res_pAN_list]

    # Summary : h2, h2s, h2_enrich, h2s_enrich
    for term in ["h2", "h2s"]:
        if term == "h2":
            v_coef = np.array(
                [dic_coef[x] for x in res_AN_list + res_pAN_list], dtype=np.float32
            )
            v_coef_jn = np.array(
                [dic_coef_jn[x] for x in res_AN_list + res_pAN_list], dtype=np.float32
            )
            mat_cov = df_coef_cov.loc[
                res_AN_list + res_pAN_list, res_AN_list + res_pAN_list
            ].values
            mat_coef_block = df_coef_block[res_AN_list + res_pAN_list].values
            dic_v = {
                x: np.concatenate([dic_AN_v[x], dic_AN_v_p[x]]) for x in res_AN_list
            }
            dic_v_ref = {
                x: np.concatenate([dic_AN_v_ref[x], dic_AN_v_p_ref[x]])
                for x in res_AN_list
            }
        if term == "h2s":
            v_coef = np.array([dic_coef[x] for x in res_AN_list], dtype=np.float32)
            v_coef_jn = np.array(
                [dic_coef_jn[x] for x in res_AN_list], dtype=np.float32
            )
            mat_cov = df_coef_cov.loc[res_AN_list, res_AN_list].values
            mat_coef_block = df_coef_block[res_AN_list].values
            dic_v = dic_AN_v
            dic_v_ref = dic_AN_v_ref

        for i_AN, AN in enumerate(res_AN_list):
            if dic_AN_type[AN] != "binary":
                continue

            n_snp_AN = dic_AN_n_snp[AN]
            n_snp_ref = dic_AN_n_snp_ref[AN]
            n_snp_dif = n_snp_ref - n_snp_AN

            # h2, h2_se
            df_sum_tau.loc[AN, "%s" % term] = (dic_v[AN] * v_coef_jn).sum()
            df_sum_tau.loc[AN, "%s_se" % term] = np.sqrt(
                dic_v[AN].dot(mat_cov).dot(dic_v[AN])
            )
            dic_jn[term][i_AN] = (dic_v[AN] * v_coef).sum()  # JN statistics for h2
            for i in range(n_jn_block):
                dic_jn["%s.jn" % term][i, i_AN] = (
                    dic_v[AN] * mat_coef_block[i, :]
                ).sum()

            # h2_enrich and h2_enrich_se via JN
            if n_snp_dif < n_snp_AN * 0.1:
                continue
            h2_ps = (dic_v[AN] * v_coef).sum() / n_snp_AN
            h2_ps_ref = (dic_v_ref[AN] * v_coef).sum() / n_snp_ref
            v_esti = [h2_ps / h2_ps_ref]
            dic_jn["%s_enrich" % term][i_AN] = (
                h2_ps / h2_ps_ref
            )  # JN statistics for h2_enrich
            mat_esti_jn = []
            for i in range(n_jn_block):
                h2_ps = (dic_v[AN] * mat_coef_block[i, :]).sum() / n_snp_AN
                h2_ps_ref = (dic_v_ref[AN] * mat_coef_block[i, :]).sum() / n_snp_ref
                mat_esti_jn.append(h2_ps / h2_ps_ref)
                dic_jn["%s_enrich.jn" % term][i, i_AN] = (
                    h2_ps / h2_ps_ref
                )  # JN statistics for h2_enrich
            v_mean_jn, mat_cov_jn = bjn(v_esti, mat_esti_jn, dic_res["v_h"])
            df_sum_tau.loc[AN, "%s_enrich" % term] = v_mean_jn[0]
            df_sum_tau.loc[AN, "%s_enrich_se" % term] = np.sqrt(mat_cov_jn[0, 0])

            # h2_enrich_p
            temp_v = (
                dic_v[AN] * (1 / n_snp_AN + 1 / n_snp_dif) - dic_v_ref[AN] / n_snp_dif
            )
            dif_ = (temp_v * v_coef_jn).sum()
            se_ = np.sqrt(temp_v.dot(mat_cov).dot(temp_v))
            df_sum_tau.loc[AN, "%s_enrich_p" % term] = gdreg.util.zsc2pval(
                dif_ / se_, option="two-sided"
            )

    # Summary : cov, cov_se, cor, cor_se
    v_coef = np.array([dic_coef[x] for x in res_pAN_list], dtype=np.float32)
    v_coef_jn = np.array([dic_coef_jn[x] for x in res_pAN_list], dtype=np.float32)
    mat_cov = df_coef_cov.loc[res_pAN_list, res_pAN_list].values
    mat_coef_block = df_coef_block[res_pAN_list].values
    #     for pAN in res_pAN_list:
    for i_pAN, pAN in enumerate(res_pAN_list):
        # cov, cov_se
        df_sum_rho.loc[pAN, "cov"] = (dic_pAN_v[pAN] * v_coef_jn).sum()
        df_sum_rho.loc[pAN, "cov_se"] = np.sqrt(
            dic_pAN_v[pAN].dot(mat_cov).dot(dic_pAN_v[pAN])
        )
        dic_jn["cov"][i_pAN] = (dic_pAN_v[pAN] * v_coef).sum()  # JN statistics for cov
        for i in range(n_jn_block):
            dic_jn["cov.jn"][i, i_pAN] = (dic_pAN_v[pAN] * mat_coef_block[i, :]).sum()

        # cor, cor_se via JN
        v_esti = [(dic_pAN_v[pAN] * v_coef).sum() / dic_pAN_var[pAN]]
        dic_jn["cor"][i_pAN] = v_esti[0]  # JN statistics for cor
        mat_esti_jn = []
        for i in range(n_jn_block):
            v_coef_block = mat_coef_block[i, :]
            mat_esti_jn.append(
                (dic_pAN_v[pAN] * v_coef_block).sum() / dic_pAN_var_block[pAN][i]
            )
            dic_jn["cor.jn"][i, i_pAN] = mat_esti_jn[-1]  # JN statistics for cor
        v_mean_jn, mat_cov_jn = bjn(v_esti, mat_esti_jn, dic_res["v_h"])
        df_sum_rho.loc[pAN, "cor"] = v_mean_jn[0]
        df_sum_rho.loc[pAN, "cor_se"] = np.sqrt(mat_cov_jn[0, 0])

        # ecov, ecov_se; ecor, ecor_se via JN
        temp_v = dic_pAN_v[pAN].copy()
        for prox in res_prox_list:
            temp_r = dic_pAN_v[prox][i_pAN] / df_sum_rho.loc[prox, "n_pair"]
            temp_v -= dic_pAN_v[prox] * temp_r  # TODO: support overlapping prox

        df_sum_rho.loc[pAN, "ecov"] = (temp_v * v_coef_jn).sum()
        df_sum_rho.loc[pAN, "ecov_se"] = np.sqrt(temp_v.dot(mat_cov).dot(temp_v))

        dic_jn["ecov"][i_pAN] = (temp_v * v_coef).sum()  # JN statistics for ecov
        for i in range(n_jn_block):
            dic_jn["ecov.jn"][i, i_pAN] = (temp_v * mat_coef_block[i, :]).sum()

        v_esti = [(temp_v * v_coef).sum() / dic_pAN_var[pAN]]
        dic_jn["ecor"][i_pAN] = v_esti[0]  # JN statistics for ecor
        mat_esti_jn = []
        for i in range(n_jn_block):
            v_coef_block = mat_coef_block[i, :]
            mat_esti_jn.append(
                (temp_v * v_coef_block).sum() / dic_pAN_var_block[pAN][i]
            )
            dic_jn["ecor.jn"][i, i_pAN] = mat_esti_jn[-1]  # JN statistics for ecor
        v_mean_jn, mat_cov_jn = bjn(v_esti, mat_esti_jn, dic_res["v_h"])
        df_sum_rho.loc[pAN, "ecor"] = v_mean_jn[0]
        df_sum_rho.loc[pAN, "ecor_se"] = np.sqrt(mat_cov_jn[0, 0])

    # Add p-value via z-score
    for term in ["tau", "h2", "h2s"]:
        temp_z = df_sum_tau[term].values / df_sum_tau["%s_se" % term].values
        df_sum_tau["%s_p" % term] = gdreg.util.zsc2pval(temp_z, option="two-sided")
    for term in ["rho", "cov", "ecov"]:
        temp_z = df_sum_rho[term].values / df_sum_rho["%s_se" % term].values
        df_sum_rho["%s_p" % term] = gdreg.util.zsc2pval(temp_z, option="two-sided")

    # Add h2_shrink, p-value via testing whether h2 - h2s = 0
    v_mean_jn, mat_cov_jn = bjn(
        dic_jn["h2"] / dic_jn["h2s"], dic_jn["h2.jn"] / dic_jn["h2s.jn"], dic_jn["v_h"]
    )
    df_sum_tau["h2_shrink"] = v_mean_jn
    df_sum_tau["h2_shrink_se"] = np.sqrt(np.diag(mat_cov_jn))
    v_mean_jn, mat_cov_jn = bjn(
        dic_jn["h2"] - dic_jn["h2s"], dic_jn["h2.jn"] - dic_jn["h2s.jn"], dic_jn["v_h"]
    )
    temp_z = v_mean_jn / np.sqrt(np.diag(mat_cov_jn))
    df_sum_tau["h2_shrink_p"] = gdreg.util.zsc2pval(temp_z, option="two-sided")

    # Add p-value for cor (same as cov) and ecor (same as ecov)
    df_sum_rho["cor_p"] = df_sum_rho["cov_p"]
    df_sum_rho["ecor_p"] = df_sum_rho["ecov_p"]

    if verbose:
        print(
            verbose_prefix + "    Completed, time=%0.1fs" % (time.time() - start_time)
        )

    return {
        "tau": df_sum_tau,
        "rho": df_sum_rho,
        "prox_list": res_prox_list,
        "dic_jn": dic_jn,
    }


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
    null_model=None,
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
    if "LD:AN:all" in df_reg:
        temp_list = ["LD:AN:all"]
    else:
        temp_list = [
            x for x in df_reg if x.startswith("LD:AN:all_") | x.startswith("LD:AN:mbin")
        ]
    v_ld = df_reg[temp_list].sum(axis=1).values.clip(min=0.1)
    if verbose:
        print(
            verbose_prefix
            + "    Use following annots for LD score: %s" % ",".join(temp_list)
        )
    # 2. Zsq variance weights : (clipped at max=10)
    #    - Z_i^2 : 1 / [ 2 (N l_i / M + 1) ^ 2 ]
    #    - Z_i Z_j : 1 / [ (N l_i / M + 1) (N l_j / M + 1) + (N l_ij / M + 1) ^ 2 ]
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

    # Regression
    mat_X = df_reg[reg_list].values.astype(np.float32)
    mat_X[:, :-1] *= n_sample_zsq
    v_y = df_reg["ZSQ"].values.astype(np.float32)

    mat_X = (mat_X.T * v_w).T
    v_y = v_y * v_w

    #     coef, coef_mean, coef_cov = reg_bjn(v_y, mat_X, dic_block)
    #     dic_res_reg = {
    #         "term": reg_list,
    #         "coef": coef,
    #         "coef_jn": coef_mean,
    #         "coef_jn_cov": coef_cov,
    #     }

    dic_jn = reg_bjn(v_y, mat_X, dic_block)
    dic_res_reg = {
        "term": reg_list,
        "coef": dic_jn["coef"],
        "coef_jn": dic_jn["coef_mean"],
        "coef_jn_cov": dic_jn["coef_cov"],
        "coef_block": dic_jn["coef_block"],
        "v_h": dic_jn["v_h"],
    }

    # Model evaluation
    dic_eval = {}
    if null_model is not None:
        if verbose:
            print(verbose_prefix + "    Eval against null: %s" % ", ".join(null_model))

        # Observed statistics: loglss, mse, mae
        mat_X = df_reg[reg_list].values.astype(np.float32)
        mat_X[:, :-1] *= n_sample_zsq
        v_y = df_reg["ZSQ"].values.astype(np.float32)
        if verbose:
            print(
                verbose_prefix
                + "    Chi2: %d zeros, imputed as %0.2e"
                % ((v_y == 0).sum(), v_y[v_y > 0].min())
            )
        v_y[v_y == 0] = v_y[v_y > 0].min()
        v_y_hat = np.dot(mat_X, dic_res_reg["coef"])
        temp_v = sp.stats.gamma.logpdf(v_y, a=0.5, scale=2 * v_y_hat.clip(min=0.1))
        dic_eval["loglss"] = (temp_v / v_ld).sum()
        dic_eval["sqe"] = ((v_y - v_y_hat) ** 2 / v_ld).sum()
        dic_eval["abe"] = (np.absolute(v_y - v_y_hat) / v_ld).sum()

        # Null loglss
        n_rep = 20
        dic_eval["loglss.null"] = np.zeros(n_rep, dtype=np.float32)
        dic_eval["sqe.null"] = np.zeros(n_rep, dtype=np.float32)
        dic_eval["abe.null"] = np.zeros(n_rep, dtype=np.float32)
        for i_rep in range(n_rep):
            np.random.seed(i_rep)
            mat_X = df_reg[reg_list].values.astype(np.float32)
            mat_X[:, :-1] *= n_sample_zsq
            v_y = df_reg["ZSQ"].values.astype(np.float32)
            for i_reg, reg in enumerate(reg_list):
                if reg not in null_model:
                    mat_X[:, i_reg] = np.random.permutation(mat_X[:, i_reg])

            mat_X = (mat_X.T * v_w).T
            v_y = v_y * v_w
            dic_jn = reg_bjn(v_y, mat_X, dic_block)

            mat_X = df_reg[reg_list].values.astype(np.float32)
            mat_X[:, :-1] *= n_sample_zsq
            v_y = df_reg["ZSQ"].values.astype(np.float32)
            v_y[v_y == 0] = v_y[v_y > 0].min()
            v_y_hat = np.dot(mat_X, dic_jn["coef"])
            temp_v = sp.stats.gamma.logpdf(v_y, a=0.5, scale=2 * v_y_hat.clip(min=0.1))
            dic_eval["loglss.null"][i_rep] = (temp_v / v_ld).sum()
            dic_eval["sqe.null"][i_rep] = ((v_y - v_y_hat) ** 2 / v_ld).sum()
            dic_eval["abe.null"][i_rep] = (np.absolute(v_y - v_y_hat) / v_ld).sum()

        for term in ["loglss", "sqe", "abe"]:
            dic_eval["%s.dif" % term] = (
                dic_eval[term] - dic_eval["%s.null" % term].mean()
            )
            dic_eval["%s.null_se" % term] = dic_eval["%s.null" % term].std()
            dic_eval["%s.dif_z" % term] = (
                dic_eval["%s.dif" % term] / dic_eval["%s.null_se" % term]
            )
            dic_eval["%s.dif_p" % term] = gdreg.util.zsc2pval(
                dic_eval["%s.dif_z" % term], option="two-sided"
            )

            if verbose:
                print(
                    verbose_prefix
                    + "    %s=%0.3e, dif=%0.2f, null_se=%0.2f, z=%0.2f, p=%0.2e"
                    % (
                        term,
                        dic_eval[term],
                        dic_eval["%s.dif" % term],
                        dic_eval["%s.null_se" % term],
                        dic_eval["%s.dif_z" % term],
                        dic_eval["%s.dif_p" % term],
                    )
                )
        dic_res_reg["eval"] = dic_eval

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
    coef : np.ndarray(dtype=np.float32)
        Estimates using full data of shape (n_regressor,).
    coef_mean : np.ndarray(dtype=np.float32)
        JN bias-corrected estimates of shape (n_regressor,).
    coef_cov : np.ndarray(dtype=np.float32)
        JN covariance of shape (n_regressor, n_regressor).
    coef_block : np.ndarray(dtype=np.float32)
        JN estimates of shape (n_block, n_regressor)
    v_h : np.ndarray(dtype=np.float32)
        n_sample / v_sample_block of shape (n_block,)
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
        print("# Call: gdreg.regress.reg_bjn")
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
    coef_mean, coef_cov = bjn(coef, coef_block, v_h)

    dic_res = {
        "coef": coef,
        "coef_mean": coef_mean,
        "coef_cov": coef_cov,
        "coef_block": coef_block,
        "v_h": v_h,
    }

    return dic_res


#     coef_mean = (-coef_block + coef).sum(axis=0) + (coef_block.T / v_h).sum(axis=1)

#     mat_tau = np.zeros([n_block, n_regressor], dtype=np.float32)
#     for i in np.arange(n_block):
#         mat_tau[i, :] = v_h[i] * coef - (v_h[i] - 1) * coef_block[i, :]

#     coef_cov = np.zeros([n_regressor, n_regressor], dtype=np.float32)
#     for i in np.arange(n_block):
#         temp_v = mat_tau[i, :] - coef_mean
#         coef_cov += np.outer(temp_v, temp_v) / n_block / (v_h[i] - 1)

#     return coef, coef_mean, coef_cov


def bjn(v_esti, mat_esti_jn, v_h):
    """Block jackknife.

    Parameters
    ----------
    v_esti : np.ndarray(dtype=np.float32)
        Estimates using all samples of shape (n_param, )
    mat_esti_jn : np.ndarray(dtype=np.float32)
        JN estimates of shape (n_block, n_param)
    v_h : np.ndarray(dtype=np.float32)
        n_sample / v_sample_block of shape (n_block, )

    Returns
    -------
    v_mean_jn : np.ndarray(dtype=np.float32)
        Jackknife bias-corrected estimates of shape (n_regressor,).
    mat_cov_jn : np.ndarray(dtype=np.float32)
        Jackknife covariance of shape (n_regressor, n_regressor).
    """

    v_esti = np.array(v_esti, dtype=np.float32)
    mat_esti_jn = np.array(mat_esti_jn, dtype=np.float32)
    if len(mat_esti_jn.shape) == 1:
        mat_esti_jn = mat_esti_jn.reshape([-1, 1])

    n_block, n_param = mat_esti_jn.shape
    v_mean_jn = (-mat_esti_jn + v_esti).sum(axis=0) + (mat_esti_jn.T / v_h).sum(axis=1)

    mat_tau = np.zeros([n_block, n_param], dtype=np.float32)
    for i in np.arange(n_block):
        mat_tau[i, :] = v_h[i] * v_esti - (v_h[i] - 1) * mat_esti_jn[i, :]

    mat_cov_jn = np.zeros([n_param, n_param], dtype=np.float32)
    for i in np.arange(n_block):
        temp_v = mat_tau[i, :] - v_mean_jn
        mat_cov_jn += np.outer(temp_v, temp_v) / n_block / (v_h[i] - 1)

    return v_mean_jn, mat_cov_jn
