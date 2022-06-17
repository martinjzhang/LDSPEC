import numpy as np
import scipy as sp
import pandas as pd
import time
import gdreg
import warnings


def simulate_snp_effect(
    df_annot,
    dic_coef,
    pannot_list=[],
    pannot_hr_list=[],
    h2g=0.5,
    alpha=-0.38,
    p_causal=0.2,
    block_size=2000,
    sym_non_pAN="non-pAN",
    random_seed=0,
    verbose=False,
):

    """
    Simulate casual SNP effects.

    Parameters
    ----------
    df_annot : pd.DataFrame
        Single-SNP annotation. Should contain all SNPs and columns
        ['CHR', 'SNP', 'BP', 'AN:name', ...]. Should also contain
        'MAF' column if alpha!=-1.
    dic_coef : dic
        GDREG model coefficients `tau` and `rho`.

        - dic[AN] : `tau` for AN `c` defined as
            `Var(\beta_i) = \sum_c a_ci tau_c`
        - dic[pAN|AN] : `rho` for pAN `l` defined as (note we simulate cor instead of cov here)
            `Cor(\beta_i, \beta_j) = \sum_{k,c} G_{k,ij} * 0.5 * [a_ci + a_cj] \rho_l`

    df_pannot_list : list of pd.DataFrame
        Each element corresponds to SNP-pair annotation. Must contain
        ['CHR', 'SNP', 'BP', 'pAN:name'] columns.
    df_pannot_hr_list : list of pd.DataFrame
        Each element corresponds to a high-res SNP-pair annotation. Must contain
        ['CHR', 'SNP', 'BP', 'pCHR', 'pSNP', 'pBP', 'pAN:name'] columns.
    h2g : float, default=0.5
        Overall heritability, equal to `sum(effect**2)`.
    alpha : float, default=-0.38
        Parameter for maf-dependent architecture,
        `Var(beta_j) \propto [maf (1-maf)]^(1+\alpha)`,
        where `\beta_j` is the standardized SNP effect size. `alpha=-1` means there
        is no maf-dependency. Schoech et al. NC 2019 suggestsed alpha=-0.38.
    p_causal : float, default=0.2
        Proportion of causal SNPs.
    block_size : int, default=2000
        Maximum number of SNPs to simulate at a time.
    sym_non_pAN : str, default='non-pAN'
        Symbol for SNPs not in the SNP-pair annotation.
    random_seed : int, default=0
        Random seed.
    verbose : bool, default=False
        If to output messages.

    Returns
    -------
    df_effect: pd.DataFrame
        Simulated SNP effects, with columns ['CHR', 'SNP', 'MAF', 'EFF'].

    """

    np.random.seed(random_seed)
    start_time = time.time()

    # Get info and add 0's for GDREG coefficients not in dic_coef
    CHR_list = sorted(set(df_annot["CHR"]))
    AN_list = [x for x in df_annot if x.startswith("AN:")]
    pAN_list = [[y for y in x if y.startswith("pAN")][0] for x in pannot_list]
    pAN_hr_list = [[y for y in x if y.startswith("pAN")][0] for x in pannot_hr_list]

    for AN in [x for x in AN_list if x not in dic_coef]:
        dic_coef[AN] = 0
    for pAN in [x for x in pAN_list if x not in dic_coef]:
        dic_coef["%s|%s" % (pAN, AN)] = 0

    # df_effect
    df_effect = df_annot.copy()
    df_effect["VAR"] = 0
    df_effect["EFF"] = 0
    df_effect.drop_duplicates("SNP", inplace=True)
    df_effect.index = df_effect["SNP"]
    df_effect.sort_values(by=["CHR", "BP"], inplace=True)
    n_snp = df_effect.shape[0]

    # Read pannot file and pannot_hr file
    dic_pannot = {}
    for temp_df in pannot_list:
        pAN = [x for x in temp_df if x.startswith("pAN:")][0]
        dic_pannot[pAN] = {x: y for x, y in zip(temp_df["SNP"], temp_df[pAN])}
    dic_pannot_hr = {}
    for temp_df in pannot_hr_list:
        pAN = [x for x in temp_df if x.startswith("pAN:")][0]
        dic_pannot_hr[pAN] = [(x, y) for x, y in zip(temp_df["SNP"], temp_df["pSNP"])]

    if verbose:
        print("# Call: gdreg.simulate.simulate_snp_effect")
        print(
            "    %d SNPs, h2g=%0.2f, alpha=%0.2f, p_causal=%0.2f"
            % (n_snp, h2g, alpha, p_causal)
        )
        print("    Single-SNP annots : %s" % ", ".join(AN_list))
        print(
            "    SNP-pair annots : %s"
            % ", ".join(pAN_list + ["%s (hr)" % x for x in pAN_hr_list])
        )
        print(
            "    Coefficients : ",
            {x: dic_coef[x] for x in dic_coef if dic_coef[x] != 0},
        )

    # Per-SNP variance + genetic architecture
    v_var = np.zeros(n_snp)
    for AN in AN_list:
        v_var = v_var + df_effect[AN] * dic_coef[AN]

    if p_causal != 1:
        # Sparsity
        ind_select = np.random.binomial(1, 1 - p_causal, size=n_snp) == 1
        v_var[ind_select] = 0

    if alpha != -1:
        # MAF-dependent architecture
        v_maf = df_effect["MAF"].values
        v_var = v_var * (v_maf * (1 - v_maf)) ** (1 + alpha)

    v_var = v_var * h2g / v_var.sum()
    df_effect["VAR"] = v_var

    # Simulate SNP causal effects one block at a time
    for CHR in CHR_list:
        # Block size
        v_snp_chr = df_effect.loc[df_effect["CHR"] == CHR, "SNP"].values
        n_snp_chr = v_snp_chr.shape[0]
        n_block_chr = int(np.ceil(n_snp_chr / block_size))
        block_size_chr = int(np.ceil(n_snp_chr / n_block_chr))
        for i in range(n_block_chr):
            v_snp_block = v_snp_chr[i * block_size_chr : (i + 1) * block_size_chr]
            n_snp_block = v_snp_block.shape[0]
            mat_annot = df_effect.loc[v_snp_block, AN_list].values.copy()

            # Correlation from SNP-pair annotations
            mat_cov = np.eye(n_snp_block)
            for pAN in pAN_list:
                v_pAN = [
                    dic_pannot[pAN][x] if x in dic_pannot[pAN] else sym_non_pAN
                    for x in v_snp_block
                ]
                mat_S = gdreg.util.pannot_to_csr(v_pAN)
                mat_G = mat_S.dot(mat_S.T).toarray()
                np.fill_diagonal(mat_G, 0)
                mat_cov += mat_G * dic_coef[pAN]

            for pAN in pAN_hr_list:
                snp_set = set(v_snp_block)
                snp_pair_list = [
                    x
                    for x in dic_pannot_hr[pAN]
                    if (x[0] in snp_set) & (x[1] in snp_set)
                ]
                mat_G = gdreg.util.pannot_hr_to_csr(
                    v_snp_block, snp_pair_list
                ).toarray()
                mat_cov += mat_G * dic_coef[pAN]

            # Scale by variance
            v_var_block = df_effect.loc[v_snp_block, "VAR"].values
            v_sd_block = np.sqrt(v_var_block)
            mat_cov = (mat_cov * v_sd_block).T * v_sd_block

            # Sample effects
            df_effect.loc[v_snp_block, "EFF"] = gdreg.util.sample_mvn(
                mat_cov, random_seed=random_seed + i
            )

    df_effect = df_effect[["CHR", "SNP", "BP", "EFF"]].copy()

    if verbose:
        print(
            "    Completed, %d SNPs simulated, time=%0.1fs"
            % (df_effect.shape[0], time.time() - start_time)
        )
    return df_effect


def summarize_snp_effect(
    df_effect,
    df_annot,
    pannot_list=[],
    pannot_hr_list=[],
    block_size=2000,
    sym_non_pAN="non-pAN",
    verbose=False,
):

    """
    Summerize the simulated SNP effects

    Parameters
    ----------
    df_effect : pd.DataFrame
        Simulated SNP effects. Must contain ['CHR', 'SNP', 'BP', 'EFF'] columns.
    df_annot : pd.DataFrame
        Single-SNP annotation. Should contain all SNPs and columns
        ['CHR', 'SNP', 'BP', 'AN:AN1', ...].
    df_pannot_list : list of pd.DataFrame
        Each element corresponds to SNP-pair annotation. Must contain
        ['CHR', 'SNP', 'BP', 'pAN:pAN1'] columns.
    df_pannot_hr_list : list of pd.DataFrame
        Each element corresponds to a high-res SNP-pair annotation. Must contain
        ['CHR', 'SNP', 'BP', 'pCHR', 'pSNP', 'pBP', 'pAN:pAN1'] columns.
    block_size : int, default=2000
        Maximum number of SNPs to simulate at a time.
    sym_non_pAN : str, default='non-pAN'
        Symbol for SNPs not in the SNP-pair annotation.
    verbose : bool, default=False
        If to output messages.

    Returns
    -------
    df_summary : pd.DataFrame
        Summary of SNP effects by annotations. One single-SNP annotation per row.
        Columns are ['n_snp', 'tau', 'h2_ps', 'h2', 'rho|pAN1', 'r2_ps|pAN1', ...]
    """

    start_time = time.time()

    # Get info
    CHR_list = sorted(set(df_annot["CHR"]))
    AN_list = [x for x in df_annot if x.startswith("AN:")]
    pAN_list = [[y for y in x if y.startswith("pAN")][0] for x in pannot_list]
    pAN_hr_list = [[y for y in x if y.startswith("pAN")][0] for x in pannot_hr_list]

    # df_snp : copy df_effect and attach df_annot
    df_snp = df_effect.copy()
    df_snp.index = df_snp["SNP"]
    df_snp.drop_duplicates("SNP", inplace=True)
    df_snp.sort_values(by=["CHR", "BP"], inplace=True)

    temp_df = df_annot[AN_list].copy()
    temp_df.index = df_annot["SNP"]
    df_snp = df_snp.join(temp_df)
    n_snp = df_snp.shape[0]

    # Read pannot file and pannot_hr file
    dic_pannot = {}
    for temp_df in pannot_list:
        pAN = [x for x in temp_df if x.startswith("pAN:")][0]
        dic_pannot[pAN] = {x: y for x, y in zip(temp_df["SNP"], temp_df[pAN])}
    dic_pannot_hr = {}
    for temp_df in pannot_hr_list:
        pAN = [x for x in temp_df if x.startswith("pAN:")][0]
        dic_pannot_hr[pAN] = [(x, y) for x, y in zip(temp_df["SNP"], temp_df["pSNP"])]

    dic_mat_G = {}
    for pAN in pAN_list:
        v_pAN = [
            dic_pannot[pAN][x] if x in dic_pannot[pAN] else sym_non_pAN
            for x in df_snp["SNP"]
        ]
        mat_G = gdreg.util.pannot_to_csr(
            v_pAN, sym_non_pAN=sym_non_pAN, flag_matS=False
        )
        dic_mat_G[pAN] = mat_G.copy()

    snp_set = set(df_snp["SNP"].values)
    for pAN in pAN_hr_list:
        snp_pair_list = [
            x for x in dic_pannot_hr[pAN] if (x[0] in snp_set) & (x[1] in snp_set)
        ]
        mat_G = gdreg.util.pannot_hr_to_csr(df_snp["SNP"].values, snp_pair_list)
        dic_mat_G[pAN] = mat_G.copy()

    if verbose:
        print("# Call: gdreg.simulate.summarize_snp_effect")
        print("    %d SNPs" % n_snp)
        print("    Single-SNP annots : %s" % ", ".join(AN_list))
        print(
            "    SNP-pair annots : %s"
            % ", ".join(pAN_list + ["%s (hr)" % x for x in pAN_hr_list])
        )

    # Summary
    v_eff = (df_snp["EFF"].values != 0) * 1
    df_sum_tau = pd.DataFrame(
        index=AN_list,
        data={
            "annot": AN_list,
            "size": [(df_snp[x] == 1).sum() for x in AN_list],
            "p_causal": [v_eff[df_snp[x] == 1].mean() for x in AN_list],
            "coef": np.nan,
            "aggeff": np.nan,
        },
    )

    df_sum_rho = pd.DataFrame(
        index=pAN_list + pAN_hr_list,
        data={
            "annot": pAN_list + pAN_hr_list,
            "size": [dic_mat_G[x].sum() for x in pAN_list + pAN_hr_list],
            "p_causal": [
                dic_mat_G[x].dot(v_eff).T.dot(v_eff) / dic_mat_G[x].sum()
                for x in pAN_list + pAN_hr_list
            ],
            "coef": np.nan,
            "aggeff": np.nan,
        },
    )

    # Summary : variance terms
    v_y = df_snp["EFF"].values ** 2
    mat_X = df_snp[AN_list].values
    df_sum_tau["coef"] = gdreg.util.reg(v_y, mat_X)

    for AN in AN_list:
        if len(set(df_snp[AN])) > 2:
            continue
        # aggeff (h2)
        temp_v = df_snp.loc[df_snp[AN] == 1, AN_list].sum(axis=0).values
        df_sum_tau.loc[AN, "aggeff"] = (temp_v * df_sum_tau["coef"]).sum()

    # Summary : covariance terms
    df_reg = None
    for CHR in CHR_list:
        v_snp_chr = df_snp.loc[df_effect["CHR"] == CHR, "SNP"].values
        n_snp_chr = v_snp_chr.shape[0]
        n_block_chr = int(np.ceil(n_snp_chr / block_size))
        block_size_chr = int(np.ceil(n_snp_chr / n_block_chr))
        for i in range(n_block_chr):
            v_snp_block = v_snp_chr[i * block_size_chr : (i + 1) * block_size_chr]
            n_snp_block = v_snp_block.shape[0]
            mat_annot = df_snp.loc[v_snp_block, AN_list].values.copy()
            temp_dic_reg = {}

            # beta_i * beta_j
            v_beta = df_snp.loc[v_snp_block, "EFF"].values
            temp_dic_reg["beta_i_beta_j"] = np.outer(v_beta, v_beta)[
                np.triu_indices(n_snp_block, k=1)
            ]

            # Regressors from .pannot
            for pAN in pAN_list:
                v_pAN = [
                    dic_pannot[pAN][x] if x in dic_pannot[pAN] else sym_non_pAN
                    for x in v_snp_block
                ]
                mat_S = gdreg.util.pannot_to_csr(v_pAN)
                mat_G = mat_S.dot(mat_S.T).toarray()
                np.fill_diagonal(mat_G, 0)
                temp_dic_reg[pAN] = mat_G[np.triu_indices(n_snp_block, k=1)] * 1

            for pAN in pAN_hr_list:
                snp_set = set(v_snp_block)
                snp_pair_list = [
                    x
                    for x in dic_pannot_hr[pAN]
                    if (x[0] in snp_set) & (x[1] in snp_set)
                ]
                mat_G = gdreg.util.pannot_hr_to_csr(
                    v_snp_block, snp_pair_list
                ).toarray()
                temp_dic_reg[pAN] = mat_G[np.triu_indices(n_snp_block, k=1)] * 1

            temp_df = pd.DataFrame(data=temp_dic_reg)
            if df_reg is None:
                df_reg = temp_df.copy()
            else:
                df_reg = pd.concat([df_reg, temp_df], axis=0)

    reg_list = [x for x in df_reg if x != "beta_i_beta_j"]
    df_sum_rho["coef"] = gdreg.util.reg(df_reg["beta_i_beta_j"], df_reg[reg_list])
    for pAN in pAN_list + pAN_hr_list:
        # aggeff (cov)
        n_snp_pair = df_sum_rho.loc[pAN, "size"]
        temp_v = np.array(
            [(dic_mat_G[pAN] * dic_mat_G[x]).sum() for x in pAN_list + pAN_hr_list]
        )
        df_sum_rho.loc[pAN, "aggeff"] = (temp_v * df_sum_rho["coef"]).sum()

    df_summary = pd.concat([df_sum_tau, df_sum_rho], axis=0)
    if verbose:
        print(df_summary)
        print("    Completed, time=%0.1fs" % (time.time() - start_time))
    return df_summary


def simulate_phen(
    df_effect,
    dic_data,
    block_size=500,
    trait_name="TRAIT",
    random_seed=0,
    verbose=False,
):
    """Simulate phenotype

    Parameters
    ----------
    df_effect : pd.DataFrame
        Simulated SNP effects. Must contain ['CHR', 'SNP', 'BP', 'EFF'] columns.
    dic_data : dict
        Genotype data reader, organized, for each CHR, as

        - dic_data[CHR]['pgen'] : .pgen file path
        - dic_data[CHR]['pvar'] : .pvar pd.DataFrame
        - dic_data[CHR]['psam'] : .psam pd.DataFrame

    block_size : int, default=500
        Maximum number of SNPs to simulate at a time.
    trait_name : str, default='TRAIT'
        Trait name
    random_seed : int, default=0
        Random seed.
    verbose : bool, default=False
        If to output messages.

    Returns
    -------
    df_phen : pd.DataFrame
        Simulated phenotypes, with columns ['FID', 'IID', 'TRAIT'].
    """

    start_time = time.time()
    np.random.seed(random_seed)

    # dic_eff
    dic_eff = {x: y for x, y in zip(df_effect["SNP"], df_effect["EFF"])}

    # dic_data : CHR_list, snp_set, sample_set, sample_set_common
    CHR_list = sorted(set(dic_data))
    snp_set = set()
    sample_set = set()
    sample_set_common = None
    for CHR in CHR_list:
        snp_set.update(dic_data[CHR]["pvar"]["SNP"])
        fid_list = dic_data[CHR]["psam"]["FID"]
        iid_list = dic_data[CHR]["psam"]["IID"]
        sample_list = ["%s:%s" % (x, y) for x, y in zip(fid_list, iid_list)]
        sample_set.update(sample_list)
        if sample_set_common is None:
            sample_set_common = sample_set.copy()
        else:
            sample_set_common = sample_set_common & set(sample_list)

    # h2g and h2e
    v_eff = np.array([dic_eff[x] if x in dic_eff else 0 for x in snp_set])
    h2g = (v_eff**2).sum()
    h2e = 1 - h2g

    if verbose:
        print("# Call: gdreg.simulate.simulate_phen")
        print(
            "    SNPs: %d in df_effect, %d in dic_data, overlap=%d"
            % (df_effect.shape[0], len(snp_set), len(set(df_effect["SNP"]) & snp_set))
        )
        print(
            "    %d/%d samples present in all %d CHRs"
            % (len(sample_set_common), len(sample_set), len(CHR_list))
        )
        print("    h2g=%0.2g, h2e=%0.2g" % (h2g, h2e))
        print("    block_size=%d, random_seed=%d" % (block_size, random_seed))

    # df_phen
    sample_list = sorted(sample_set)
    df_phen = pd.DataFrame(
        index=sample_list,
        data={
            "FID": [x.split(":")[0] for x in sample_list],
            "IID": [x.split(":")[1] for x in sample_list],
        },
    )

    # Compute phen
    for CHR in CHR_list:
        n_sample_chr = dic_data[CHR]["psam"].shape[0]
        v_phen = np.zeros(n_sample_chr, dtype=np.float32)

        n_snp_chr = dic_data[CHR]["pvar"].shape[0]
        n_block_chr = int(np.ceil(n_snp_chr / block_size))
        block_size_chr = int(np.ceil(n_snp_chr / n_block_chr))

        for i in range(n_block_chr):
            ind_start = i * block_size_chr
            ind_end = min((i + 1) * block_size_chr, n_snp_chr)
            mat_X = gdreg.util.read_geno(dic_data[CHR]["pgen"], ind_start, ind_end)
            v_snp_block = dic_data[CHR]["pvar"]["SNP"].values[ind_start:ind_end]
            v_eff = np.array([dic_eff[x] if x in dic_eff else 0 for x in v_snp_block])

            mat_X = mat_X.T.astype(np.float32)
            mat_X[mat_X == -9] = 0
            v_maf = mat_X.mean(axis=0) * 0.5
            mat_X = (mat_X - 2 * v_maf) / np.sqrt(2 * v_maf * (1 - v_maf))
            v_phen += mat_X.dot(v_eff)

        fid_list = dic_data[CHR]["psam"]["FID"]
        iid_list = dic_data[CHR]["psam"]["IID"]
        sample_list_chr = ["%s:%s" % (x, y) for x, y in zip(fid_list, iid_list)]
        temp_df = pd.DataFrame(index=sample_list_chr, data={CHR: v_phen})
        df_phen = df_phen.join(temp_df)

    v_e = np.random.randn(df_phen.shape[0]).astype(np.float32) * np.sqrt(h2e)
    df_phen["TRAIT"] = df_phen[CHR_list].sum(axis=1) + v_e
    df_phen = df_phen.loc[sorted(sample_set_common), ["FID", "IID", "TRAIT"]].copy()

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
        Genotype data reader, organized, for each CHR, as

        - dic_data[CHR]['pgen'] : .pgen file path
        - dic_data[CHR]['pvar'] : .pvar pd.DataFrame
        - dic_data[CHR]['psam'] : .psam pd.DataFrame

    block_size : int, default=500
        Maximum number of SNPs to simulate at a time.
    verbose : bool, default=False
        If to output messages.

    Returns
    -------
    df_sumstats : pd.DataFrame
        Summary statistics, with columns ['SNP', 'N', 'Z', 'A1', 'A2'].
    """

    start_time = time.time()

    # dic_phen
    phen_name = df_phen.columns[2]
    dic_phen = {
        "%s:%s" % (x, y): z
        for x, y, z in zip(df_phen["FID"], df_phen["IID"], df_phen[phen_name])
    }

    # dic_data : CHR_list, snp_set, sample_set, sample_set_common
    CHR_list = sorted(set(dic_data))
    snp_set = set()
    sample_set_common = set(dic_phen.keys())
    for CHR in CHR_list:
        snp_set.update(dic_data[CHR]["pvar"]["SNP"])
        fid_list = dic_data[CHR]["psam"]["FID"]
        iid_list = dic_data[CHR]["psam"]["IID"]
        sample_list = ["%s:%s" % (x, y) for x, y in zip(fid_list, iid_list)]
        sample_set_common = sample_set_common & set(sample_list)

    if verbose:
        print("# Call: gdreg.simulate.compute_sumstats")
        print("    %d SNPs in dic_data" % len(snp_set))
        print(
            "    %d samples in df_phen, %d/%d present in all %d CHRs in dic_data"
            % (len(dic_phen), len(sample_set_common), len(dic_phen), len(CHR_list))
        )
        print("    block_size=%d" % (block_size))

    # Compute phen
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
            mat_X[mat_X == -9] = 0
            v_maf = mat_X.mean(axis=0) * 0.5
            mat_X = (mat_X - 2 * v_maf) / np.sqrt(2 * v_maf * (1 - v_maf))

            v_sumstats.extend(v_phen.T.dot(mat_X) / np.sqrt(N_chr))

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
