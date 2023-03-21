import numpy as np
import scipy as sp
import pandas as pd
import time
import gdreg
import warnings


def simulate_snp_effect(
    dic_data,
    dic_coef,
    dic_annot_path={},
    dic_pannot_path={},
    h2g=0.5,
    alpha=-0.38,
    p_causal=0.2,
    block_size=100,
    flag_bw_sparse=False,
    random_seed=0,
    verbose=False,
):

    """
    Simulate casual SNP effects.

    Parameters
    ----------
    dic_data : dict
        Genotype data reader. Must contain SNPs from all chromosomes.

        - dic_data[CHR]['pgen'] : .pgen file path
        - dic_data[CHR]['pvar'] : .pvar pd.DataFrame
        - dic_data[CHR]['psam'] : .psam pd.DataFrame
        - dic_data[CHR]['afreq'] : .psam pd.DataFrame

    dic_coef : dic
        GDREG model coefficients `tau` and `rho`. Should match annots in dic_annot_path
        and pannots in dic_pannot_path

        - dic[AN] : `tau` for AN `c` defined as
            `Var(\beta_i) = \sum_c a_ci tau_c`
        - dic[pAN] : `rho` for pAN `k` defined as
            `Cov(\beta_i, \beta_j) = \sum_{k} G^{(k)}_{ij} \rho_k`

    dic_annot_path : dic of dic of strs
        File path for single-SNP annotation. dic_annot_path[annot_name][CHR] contains the
        `.annot.gz` file path for annot_file `annot_name` and CHR `CHR`.
    dic_pannot_path : dic of dic of strs
        File path for SNP-pair annotation. dic_pannot_path[annot_name][CHR] contains the
        `.pannot_mat.npz` file path for annotation pAN and and CHR `CHR`. Dimension of the
        sparse matrix should match `dic_data[CHR][pvar]`.
    h2g : float, default=0.5
        Overall heritability, equal to `sum(effect**2)`.
    alpha : float, default=-0.38
        Parameter for maf-dependent architecture,
        `Var(beta_j) \propto [maf (1-maf)]^(1+\alpha)`,
        where `\beta_j` is the standardized SNP effect size. `alpha=-1` means there
        is no maf-dependency. Schoech et al. NC 2019 suggestsed alpha=-0.38.
    p_causal : float, default=0.2
        Proportion of causal SNPs.
    block_size : int, default=100
        Maximum number of SNPs to simulate at a time.
    flag_bw_sparse : bool, default=False
        If True, each block (with size `block_size`) is either all causal
        (with probability `p_causal`) or all non-causal

    Returns
    -------
    df_effect: pd.DataFrame
        Simulated SNP effects, with columns ['CHR', 'SNP', 'MAF', 'EFF'].

    """

    np.random.seed(random_seed)
    start_time = time.time()

    CHR_list = list(dic_data)
    AN_list = [x for x in dic_coef if x.startswith("AN:")]
    pAN_list = [x for x in dic_coef if x.startswith("pAN:")]

    temp_list, CHR0 = [], CHR_list[0]
    for annot_name in dic_annot_path:
        temp_df = gdreg.util.read_annot(dic_annot_path[annot_name][CHR0], nrows=5)
        temp_list.extend(temp_df)
    for AN in AN_list:
        if AN not in temp_list:
            raise ValueError("%s not in dic_annot_path" % AN)
    for pAN in pAN_list:
        if pAN not in dic_pannot_path:
            raise ValueError("%s not in dic_pannot_path" % pAN)

    if verbose:
        print("# Call: gdreg.simulate.simulate_snp_effect")
        print("    h2g=%0.2f, alpha=%0.2f, p_causal=%0.2f" % (h2g, alpha, p_causal))
        print("    block_size=%d, flag_bw_sparse=%s" % (block_size, flag_bw_sparse))
        print(
            "    annots: %s"
            % ", ".join(["%s (%0.2f)" % (x, dic_coef[x]) for x in AN_list])
        )
        print(
            "    pannots: %s"
            % ", ".join(["%s (%0.2f)" % (x, dic_coef[x]) for x in pAN_list])
        )

    # Simulate
    df_list = []
    for CHR in CHR_list:
        # df_effect_chr : 'MAF', 'VAR', 'EFF', annots
        df_effect_chr = dic_data[CHR]["pvar"].copy()
        n_snp_chr = df_effect_chr.shape[0]
        df_effect_chr["MAF"] = dic_data[CHR]["afreq"]["MAF"]
        df_effect_chr["VAR"] = 0
        df_effect_chr["EFF"] = 0

        for annot_name in dic_annot_path:
            temp_df = gdreg.util.read_annot(dic_annot_path[annot_name][CHR])
            for AN in AN_list:
                if AN in temp_df:
                    temp_dic = {
                        x: y for x, y in zip(temp_df["SNP"], temp_df[AN]) if y != 0
                    }
                    df_effect_chr[AN] = [
                        temp_dic[x] if x in temp_dic else 0
                        for x in df_effect_chr["SNP"]
                    ]

        for col in ["MAF", "VAR", "EFF"] + AN_list:
            df_effect_chr[col] = df_effect_chr[col].astype(np.float32)

        # Per-SNP variance : df_effect_chr["VAR"]
        for AN in AN_list:
            df_effect_chr["VAR"] += df_effect_chr[AN] * dic_coef[AN]
        if verbose:
            print(
                "    CHR%d, VAR_min=%0.2e, %d/%d negative"
                % (
                    CHR,
                    df_effect_chr["VAR"].min(),
                    (df_effect_chr["VAR"] < 0).sum(),
                    df_effect_chr.shape[0],
                )
            )
        df_effect_chr["VAR"] = df_effect_chr["VAR"].clip(lower=0)

        # Each SNP has independent probability to be causal
        if (p_causal != 1) & (flag_bw_sparse is False):
            ind_select = np.random.binomial(1, p_causal, size=n_snp_chr) == 0
            df_effect_chr.loc[ind_select, "VAR"] = 0

        if alpha != -1:
            v_maf = df_effect_chr["MAF"].values
            df_effect_chr["VAR"] *= (v_maf * (1 - v_maf)) ** (1 + alpha)

        # Read pannot file
        dic_mat_G_chr = {}
        for pAN in pAN_list:
            dic_mat_G_chr[pAN] = gdreg.util.read_pannot_mat(dic_pannot_path[pAN][CHR])

        # Simulate SNP causal effects one block at a time
        n_block_chr = int(np.ceil(n_snp_chr / block_size))
        block_size_chr = int(np.ceil(n_snp_chr / n_block_chr))
        for i in range(n_block_chr):
            if i % 100 == 0:
                print(
                    "    CHR%2d block %d/%d, time=%0.1fs"
                    % (CHR, i, n_block_chr, time.time() - start_time)
                )

            # Block-wise sparsity : each block is causal with probability p_causal
            if flag_bw_sparse & (np.random.rand(1)[0] > p_causal):
                continue

            ind_block = np.zeros(n_snp_chr, dtype=bool)
            ind_block[i * block_size_chr : (i + 1) * block_size_chr] = True
            n_snp_block = ind_block.sum()

            # Correlation from SNP-pair annotations
            mat_cov = np.eye(n_snp_block, dtype=np.float32)
            for pAN in pAN_list:
                temp_mat_G = dic_mat_G_chr[pAN][ind_block, :].toarray()[:, ind_block]
                mat_cov += temp_mat_G * dic_coef[pAN]
            mat_cov = mat_cov.clip(-1, 1)

            # Scale by variance
            v_sd_block = df_effect_chr["VAR"].values[ind_block]
            v_sd_block = np.sqrt(v_sd_block)
            mat_cov = (mat_cov * v_sd_block).T * v_sd_block

            # Sample effects
            # df_effect_chr.loc[ind_block, "EFF"] = gdreg.util.sample_mvn(
            #     mat_cov, random_seed=random_seed + i + 500 * CHR
            # )
            ind_select = np.diag(mat_cov) > 0
            temp_v = np.zeros(n_snp_block, dtype=np.float32)
            temp_v[ind_select] = gdreg.util.sample_mvn(
                mat_cov[ind_select, :][:, ind_select]
            )
            # # Sanity check: only for diagonal elements
            # ind_select = np.diag(mat_cov)>0
            # temp_v = np.zeros(n_snp_block, dtype=np.float32)
            # temp_v[ind_select] = np.sqrt(np.diag(mat_cov)[ind_select]) * np.random.randn(ind_select.shape[0])
            df_effect_chr.loc[ind_block, "EFF"] = temp_v

        df_list.append(df_effect_chr)

    df_effect = pd.concat(df_list, axis=0)
    h2g_empi = (df_effect["EFF"] ** 2).sum()
    df_effect["EFF"] = df_effect["EFF"] * np.sqrt(h2g / h2g_empi)
    df_effect = df_effect[["CHR", "SNP", "BP", "EFF"]].copy()
    df_effect["EFF"] = df_effect["EFF"].astype(np.float32)

    if verbose:
        print(
            "    h2g=%0.3f,  p_causal=%0.3f"
            % ((df_effect["EFF"] ** 2).sum(), (df_effect["EFF"] != 0).mean())
        )
        print(
            "    Completed, %d SNPs simulated, time=%0.1fs"
            % (df_effect.shape[0], time.time() - start_time)
        )
    return df_effect


def summarize_snp_effect(
    dic_data,
    dic_coef,
    df_effect,
    #     df_phen,
    df_phen=None,
    dic_annot_path={},
    dic_pannot_path={},
    block_size=1000,
    verbose=False,
):

    """
    Summerize the simulated SNP effects

    Parameters
    ----------
    dic_data : dict
        Genotype data reader. Must contain SNPs from all chromosomes.

        - dic_data[CHR]['pgen'] : .pgen file path
        - dic_data[CHR]['pvar'] : .pvar pd.DataFrame
        - dic_data[CHR]['psam'] : .psam pd.DataFrame
        - dic_data[CHR]['afreq'] : .psam pd.DataFrame

    dic_coef : dic
        GDREG model coefficients `tau` and `rho`. Should match annots in dic_annot_path
        and pannots in dic_pannot_path

        - dic[AN] : `tau` for AN `c` defined as
            `Var(\beta_i) = \sum_c a_ci tau_c`
        - dic[pAN] : `rho` for pAN `k` defined as
            `Cov(\beta_i, \beta_j) = \sum_{k} G^{(k)}_{ij} \rho_k`

    df_effect : pd.DataFrame
        Simulated SNP effects. Must contain ['CHR', 'SNP', 'BP', 'EFF'] columns.
    df_phen : pd.DataFrame, default=None
        Simulated phenotypes, with columns ['FID', 'IID', 'TRAIT'].
    dic_annot_path : dic of dic of strs
        File path for single-SNP annotation. dic_annot_path[annot_name][CHR] contains the
        `.annot.gz` file path for annot_file `annot_name` and CHR `CHR`.
    dic_pannot_path : dic of dic of strs
        File path for SNP-pair annotation. dic_pannot_path[annot_name][CHR] contains the
        `.pannot_mat.npz` file path for annotation pAN and and CHR `CHR`. Dimension of the
        sparse matrix should match `dic_data[CHR][pvar]`.
    block_size : int, default=1000
        Maximum number of SNPs to simulate at a time.

    Returns
    -------
    df_sum_tau : pd.DataFrame
        Summary of single-SNP results. One annotation per row including columns
        ['annot', 'n_snp', 'p_causal', 'tau']
    df_sum_rho : pd.DataFrame
        Summary of SNP-pair results. One annotation per row including columns
        ['pannot', 'n_pair', 'p_causal', 'rho', 'cov', 'cor']
    """

    start_time = time.time()

    # Get info
    CHR_list = list(dic_data)
    AN_list = [x for x in dic_coef if x.startswith("AN:")]
    pAN_list = [x for x in dic_coef if x.startswith("pAN:")]
    dic_eff = {x: y for x, y in zip(df_effect["SNP"], df_effect["EFF"])}
    prox_list = [x for x in pAN_list if "prox" in x]  # pAN's for priximity

    if verbose:
        print("# Call: gdreg.simulate.summarize_snp_effect")
        print(
            "    annots: %s"
            % ", ".join(["%s (%0.2f)" % (x, dic_coef[x]) for x in AN_list])
        )
        print(
            "    pannots: %s"
            % ", ".join(["%s (%0.2f)" % (x, dic_coef[x]) for x in pAN_list])
        )
        print("    prox: %s" % ", ".join(prox_list))

    # Summary
    df_sum_tau = pd.DataFrame(
        index=AN_list,
        data={
            "annot": AN_list,
            "n_snp": np.nan,
            "p_causal": np.nan,
            "tau": np.nan,
            "h2": np.nan,
            "h2_enrich": np.nan,
        },
    )

    df_sum_rho = pd.DataFrame(
        index=pAN_list,
        data={
            "pannot": pAN_list,
            "n_pair": 0,
            "p_causal": 0,
            "rho": np.nan,
            "cov": 0,
            "cor": 0,
            "ecov": 0,
            "ecor": 0,
        },
    )

    # Summary : single-SNP annotations
    df_list = []
    for CHR in CHR_list:
        df_snp_chr = dic_data[CHR]["pvar"].copy()
        n_snp_chr = df_snp_chr.shape[0]
        df_snp_chr["EFF"] = [dic_eff[x] for x in df_snp_chr["SNP"]]

        # Read annot file
        for annot_name in dic_annot_path:
            temp_df = gdreg.util.read_annot(dic_annot_path[annot_name][CHR])
            ref_col_list = [x for x in temp_df if x.endswith("_common")]
            temp_df["AN:all_common"] = (temp_df[ref_col_list].values == 1).sum(
                axis=1
            ) > 0
            ref_col_list = [x for x in temp_df if x.endswith("_lf")]
            temp_df["AN:all_lf"] = (temp_df[ref_col_list].values == 1).sum(axis=1) > 0
            for AN in AN_list + ["AN:all_lf", "AN:all_common"]:
                if AN in temp_df:
                    temp_dic = {
                        x: y for x, y in zip(temp_df["SNP"], temp_df[AN]) if y != 0
                    }
                    df_snp_chr[AN] = [
                        temp_dic[x] if x in temp_dic else 0 for x in df_snp_chr["SNP"]
                    ]

        for col in ["EFF"] + AN_list + ["AN:all_lf", "AN:all_common"]:
            df_snp_chr[col] = df_snp_chr[col].astype(np.float32)
        df_list.append(df_snp_chr.copy())

    df_snp = pd.concat(df_list, axis=0)
    v_y = df_snp["EFF"].values ** 2
    mat_X = df_snp[AN_list].values
    df_sum_tau["tau"] = gdreg.util.reg(v_y, mat_X)
    df_sum_tau["n_snp"] = [(df_snp[x] == 1).sum() for x in AN_list]
    df_sum_tau["p_causal"] = [
        (df_snp.loc[df_snp[x] == 1, "EFF"] != 0).mean() for x in AN_list
    ]

    if df_phen is not None:
        h2_common = df_phen["AN:all_common"].var() / df_phen["TRAIT"].var()
        h2_common_ps = h2_common / (df_snp["AN:all_common"] == 1).sum()
        h2_lf = df_phen["AN:all_lf"].var() / df_phen["TRAIT"].var()
        h2_lf_ps = h2_lf / (df_snp["AN:all_lf"] == 1).sum()
        for AN in AN_list:
            df_sum_tau.loc[AN, "h2"] = df_phen[AN].var() / df_phen["TRAIT"].var()
            h2ps = df_sum_tau.loc[AN, "h2"] / df_sum_tau.loc[AN, "n_snp"]
            if AN.endswith("_common"):
                df_sum_tau.loc[AN, "h2_enrich"] = h2ps / h2_common_ps
            if AN.endswith("_lf"):
                df_sum_tau.loc[AN, "h2_enrich"] = h2ps / h2_lf_ps
    else:
        df_sum_tau["h2_enrich"] = 0

    v_h2ps = np.zeros(df_snp.shape[0], dtype=np.float32)
    for AN in AN_list:
        v_h2ps += df_sum_tau.loc[AN, "tau"] * df_snp[AN].values
    v_h2ps = np.sqrt(v_h2ps.clip(min=0)).astype(np.float32)
    dic_h2ps = {x: y for x, y in zip(df_snp["SNP"], v_h2ps)}

    # Summary : SNP-pair annotations
    df_list = []
    dic_overlap = {x: {y: 0 for y in prox_list} for x in pAN_list}
    for CHR in CHR_list:
        df_snp_chr = dic_data[CHR]["pvar"].copy()
        n_snp_chr = df_snp_chr.shape[0]
        df_snp_chr["EFF"] = [dic_eff[x] for x in df_snp_chr["SNP"]]
        v_h2ps_chr = np.array(
            [dic_h2ps[x] for x in df_snp_chr["SNP"]], dtype=np.float32
        )

        # Read pannot file
        dic_mat_G_chr = {}
        for pAN in pAN_list:
            dic_mat_G_chr[pAN] = gdreg.util.read_pannot_mat(dic_pannot_path[pAN][CHR])
            df_sum_rho.loc[pAN, "n_pair"] += dic_mat_G_chr[pAN].sum()

            v_eff_b = (df_snp_chr["EFF"].values != 0) * 1
            df_sum_rho.loc[pAN, "p_causal"] += (
                dic_mat_G_chr[pAN].dot(v_eff_b).T.dot(v_eff_b)
            )

            v_eff = df_snp_chr["EFF"].values
            df_sum_rho.loc[pAN, "cov"] += dic_mat_G_chr[pAN].dot(v_eff).T.dot(v_eff)
            df_sum_rho.loc[pAN, "cor"] += (
                dic_mat_G_chr[pAN].dot(v_h2ps_chr).T.dot(v_h2ps_chr)
            )

        for pAN in pAN_list:
            for prox in prox_list:
                dic_overlap[pAN][prox] += (
                    dic_mat_G_chr[pAN].multiply(dic_mat_G_chr[prox]).sum()
                )

        n_block_chr = int(np.ceil(n_snp_chr / block_size))
        block_size_chr = int(np.ceil(n_snp_chr / n_block_chr))
        for i in range(n_block_chr):
            ind_block = np.zeros(n_snp_chr, dtype=bool)
            ind_block[i * block_size_chr : (i + 1) * block_size_chr] = True
            n_snp_block = ind_block.sum()

            temp_dic_reg = {}
            # beta_i * beta_j
            v_beta = df_snp_chr["EFF"].values[ind_block]
            temp_dic_reg["beta_ij"] = np.outer(v_beta, v_beta)[
                np.triu_indices(n_snp_block, k=1)
            ]

            # Regressors from .pannot
            for pAN in pAN_list:
                temp_mat_G = dic_mat_G_chr[pAN][ind_block, :].toarray()[:, ind_block]
                temp_dic_reg[pAN] = temp_mat_G[np.triu_indices(n_snp_block, k=1)] * 1

            temp_df = pd.DataFrame(data=temp_dic_reg)
            for col in ["beta_ij"] + pAN_list:
                temp_df[col] = temp_df[col].astype(np.float32)
            df_list.append(temp_df.copy())

    df_reg = pd.concat(df_list, axis=0)
    df_sum_rho["rho"] = gdreg.util.reg(df_reg["beta_ij"], df_reg[pAN_list])
    df_sum_rho["p_causal"] = df_sum_rho["p_causal"] / df_sum_rho["n_pair"]
    df_sum_rho["ecov"] = df_sum_rho["cov"]
    for pAN in pAN_list:
        for prox in prox_list:
            temp_r = dic_overlap[pAN][prox] / df_sum_rho.loc[prox, "n_pair"]
            df_sum_rho.loc[pAN, "ecov"] -= temp_r * df_sum_rho.loc[prox, "cov"]
    temp_v = df_sum_rho["cor"].values.copy()
    df_sum_rho["cor"] = df_sum_rho["cov"] / temp_v
    df_sum_rho["ecor"] = df_sum_rho["ecov"] / temp_v

    if verbose:
        print(df_sum_tau)
        print(df_sum_rho)
        print("    Completed, time=%0.1fs" % (time.time() - start_time))
    return df_sum_tau, df_sum_rho


def simulate_phen(
    dic_data,
    dic_coef,
    df_effect,
    dic_annot_path={},
    block_size=500,
    random_seed=0,
    verbose=False,
):
    """Simulate phenotype

    Parameters
    ----------
    dic_data : dict
        Genotype data reader. Must contain SNPs from all chromosomes.

        - dic_data[CHR]['pgen'] : .pgen file path
        - dic_data[CHR]['pvar'] : .pvar pd.DataFrame
        - dic_data[CHR]['psam'] : .psam pd.DataFrame
        - dic_data[CHR]['afreq'] : .psam pd.DataFrame

    dic_coef : dic
        GDREG model coefficients `tau` and `rho`. Should match annots in dic_annot_path
        and pannots in dic_pannot_path

        - dic[AN] : `tau` for AN `c` defined as
            `Var(\beta_i) = \sum_c a_ci tau_c`
        - dic[pAN] : `rho` for pAN `k` defined as
            `Cov(\beta_i, \beta_j) = \sum_{k} G^{(k)}_{ij} \rho_k`

    df_effect : pd.DataFrame
        Simulated SNP effects. Must contain ['CHR', 'SNP', 'BP', 'EFF'] columns.
    dic_annot_path : dic of dic of strs
        File path for single-SNP annotation. dic_annot_path[annot_name][CHR] contains the
        `.annot.gz` file path for annot_file `annot_name` and CHR `CHR`.
    block_size : int, default=500
        Maximum number of SNPs to simulate at a time.

    Returns
    -------
    df_phen : pd.DataFrame
        Simulated phenotypes, with columns ['FID', 'IID', 'TRAIT'].
    """

    start_time = time.time()
    np.random.seed(random_seed)

    # Get info
    CHR_list = list(dic_data)
    AN_list = [x for x in dic_coef if x.startswith("AN:")]
    pAN_list = [x for x in dic_coef if x.startswith("pAN:")]
    dic_eff = {x: y for x, y in zip(df_effect["SNP"], df_effect["EFF"])}

    fid_list = dic_data[CHR_list[0]]["psam"]["FID"]
    iid_list = dic_data[CHR_list[0]]["psam"]["IID"]

    snp_list = []
    for CHR in CHR_list:
        snp_list.extend(dic_data[CHR]["pvar"]["SNP"])

    h2g = (df_effect["EFF"] ** 2).sum()
    h2e = 1 - h2g

    if verbose:
        print("# Call: gdreg.simulate.simulate_phen")
        print(
            "    %d samples, %d SNPs (%d in df_effect)"
            % (len(fid_list), len(snp_list), len(set(snp_list) & set(df_effect["SNP"])))
        )
        print("    h2g=%0.2g, h2e=%0.2g" % (h2g, h2e))
        print("    block_size=%d, random_seed=%d" % (block_size, random_seed))

    # df_phen
    df_phen = pd.DataFrame(
        data={
            "FID": fid_list,
            "IID": iid_list,
        },
    )
    for AN in ["TRAIT"] + AN_list + ["AN:all_common", "AN:all_lf"]:
        df_phen[AN] = 0
        df_phen[AN] = df_phen[AN].astype(np.float32)

    for CHR in CHR_list:
        # df_snp_chr with annots
        df_snp_chr = dic_data[CHR]["pvar"].copy()
        n_snp_chr = df_snp_chr.shape[0]

        for annot_name in dic_annot_path:
            temp_df = gdreg.util.read_annot(dic_annot_path[annot_name][CHR])
            ref_col_list = [x for x in temp_df if x.endswith("_common")]
            temp_df["AN:all_common"] = (temp_df[ref_col_list].values == 1).sum(
                axis=1
            ) > 0
            ref_col_list = [x for x in temp_df if x.endswith("_lf")]
            temp_df["AN:all_lf"] = (temp_df[ref_col_list].values == 1).sum(axis=1) > 0
            for AN in AN_list + ["AN:all_common", "AN:all_lf"]:
                if AN in temp_df:
                    temp_dic = {
                        x: y for x, y in zip(temp_df["SNP"], temp_df[AN]) if y != 0
                    }
                    df_snp_chr[AN] = [
                        temp_dic[x] if x in temp_dic else 0 for x in df_snp_chr["SNP"]
                    ]

        for col in AN_list + ["AN:all_common", "AN:all_lf"]:
            df_snp_chr[col] = df_snp_chr[col].astype(np.float32)

        n_block_chr = int(np.ceil(n_snp_chr / block_size))
        block_size_chr = int(np.ceil(n_snp_chr / n_block_chr))

        for i in range(n_block_chr):
            ind_start = i * block_size_chr
            ind_end = min((i + 1) * block_size_chr, n_snp_chr)
            mat_X = gdreg.util.read_geno(dic_data[CHR]["pgen"], ind_start, ind_end)
            v_snp_block = dic_data[CHR]["pvar"]["SNP"].values[ind_start:ind_end]
            v_eff = np.array([dic_eff[x] if x in dic_eff else 0 for x in v_snp_block])

            mat_X = mat_X.T.astype(np.float32)
            mat_X[mat_X == -9] = np.nan  # Imputation by mean genotype
            v_maf = np.nanmean(mat_X, axis=0) * 0.5  # Imputation by mean genotype
            #             mat_X[mat_X == -9] = 0
            #             v_maf = mat_X.mean(axis=0) * 0.5
            mat_X = (mat_X - 2 * v_maf) / np.sqrt(2 * v_maf * (1 - v_maf))
            mat_X[np.isnan(mat_X)] = 0

            df_phen["TRAIT"] += mat_X.dot(v_eff)
            for AN in AN_list + ["AN:all_common", "AN:all_lf"]:
                ind_AN = df_snp_chr[AN][ind_start:ind_end] == 1
                if ind_AN.sum() == 0:
                    continue
                df_phen[AN] += mat_X[:, ind_AN].dot(v_eff[ind_AN])

            if i % 25 == 0:
                print(
                    "    CHR%2d block %2d/%2d, time=%0.1fs"
                    % (CHR, i, n_block_chr, time.time() - start_time)
                )

    df_phen["TRAIT"] += np.random.randn(df_phen.shape[0]) * np.sqrt(h2e)
    for col in ["TRAIT"] + AN_list + ["AN:all_common", "AN:all_lf"]:
        df_phen[col] = df_phen[col].astype(np.float32)

    if verbose:
        print("    Completed, time=%0.1fs" % (time.time() - start_time))

    return df_phen


def compute_sumstats(df_phen, dic_data, block_size=500, verbose=False):

    """Compute summary statistics

    Parameters
    ----------
    df_phen : pd.DataFrame
        Simulated phenotypes, with columns ['FID', 'IID', 'TRAIT'].
    dic_data : dict
        Genotype data reader. Must contain SNPs from all chromosomes.

        - dic_data[CHR]['pgen'] : .pgen file path
        - dic_data[CHR]['pvar'] : .pvar pd.DataFrame
        - dic_data[CHR]['psam'] : .psam pd.DataFrame
        - dic_data[CHR]['afreq'] : .psam pd.DataFrame

    block_size : int, default=500
        Maximum number of SNPs to simulate at a time.

    Returns
    -------
    df_sumstats : pd.DataFrame
        Summary statistics, with columns ['SNP', 'N', 'Z', 'A1', 'A2'].
    """

    start_time = time.time()

    # df_phen
    phen_name = df_phen.columns[2]
    dic_phen = {
        "%s:%s" % (x, y): z
        for x, y, z in zip(df_phen["FID"], df_phen["IID"], df_phen[phen_name])
    }

    # dic_data
    CHR_list = list(dic_data)
    snp_list = []
    for CHR in CHR_list:
        snp_list.extend(dic_data[CHR]["pvar"]["SNP"])

    fid_list = dic_data[CHR_list[0]]["psam"]["FID"]
    iid_list = dic_data[CHR_list[0]]["psam"]["IID"]
    sample_list = ["%s:%s" % (x, y) for x, y in zip(fid_list, iid_list)]
    sample_set_common = set(dic_phen.keys()) & set(sample_list)

    if verbose:
        print("# Call: gdreg.simulate.compute_sumstats")
        print("    dic_data: %d CHRs" % len(CHR_list))
        print("    dic_data: %d SNPs" % len(snp_list))
        print("    dic_data: %d samples" % len(sample_list))
        print(
            "    df_phen: %d/%d samples also in dic_data"
            % (len(sample_set_common), len(dic_phen))
        )
        print("    block_size=%d" % (block_size))

    # Compute sumstats
    df_sumstats = pd.DataFrame(columns=["SNP", "N", "Z", "A1", "A2"])
    for CHR in CHR_list:
        # ind_sample and v_phen
        fid_list = dic_data[CHR]["psam"]["FID"]
        iid_list = dic_data[CHR]["psam"]["IID"]
        sample_list_chr = ["%s:%s" % (x, y) for x, y in zip(fid_list, iid_list)]
        ind_sample = [x in dic_phen for x in sample_list_chr]
        v_phen = np.array(
            [dic_phen[x] if x in dic_phen else 0 for x in sample_list_chr],
            dtype=np.float32,
        )[ind_sample]
        N_chr = v_phen.shape[0]

        n_snp_chr = dic_data[CHR]["pvar"].shape[0]
        n_block_chr = int(np.ceil(n_snp_chr / block_size))
        block_size_chr = int(np.ceil(n_snp_chr / n_block_chr))

        v_sumstats = []
        for i in range(n_block_chr):
            ind_start = i * block_size_chr
            ind_end = min((i + 1) * block_size_chr, n_snp_chr)
            mat_X = gdreg.util.read_geno(dic_data[CHR]["pgen"], ind_start, ind_end)

            mat_X = mat_X[:, ind_sample].copy()
            mat_X = mat_X.T.astype(np.float32)
            mat_X[mat_X == -9] = np.nan  # Imputation by mean genotype
            v_maf = np.nanmean(mat_X, axis=0) * 0.5  # Imputation by mean genotype
            #             mat_X[mat_X == -9] = 0
            #             v_maf = mat_X.mean(axis=0) * 0.5
            mat_X = (mat_X - 2 * v_maf) / np.sqrt(2 * v_maf * (1 - v_maf))
            mat_X[np.isnan(mat_X)] = 0

            v_sumstats.extend(v_phen.T.dot(mat_X) / np.sqrt(N_chr))
            if i % 25 == 0:
                print(
                    "    CHR%2d block %2d/%2d, time=%0.1fs"
                    % (CHR, i, n_block_chr, time.time() - start_time)
                )

        temp_df = pd.DataFrame(
            data={
                "SNP": dic_data[CHR]["pvar"]["SNP"],
                "N": N_chr,
                "Z": np.array(v_sumstats, dtype=np.float32),
                "A1": dic_data[CHR]["pvar"]["ALT"],
                "A2": dic_data[CHR]["pvar"]["REF"],
            }
        )
        df_sumstats = pd.concat([df_sumstats, temp_df], axis=0)

    if verbose:
        print("    Completed, time=%0.1fs" % (time.time() - start_time))
    return df_sumstats
