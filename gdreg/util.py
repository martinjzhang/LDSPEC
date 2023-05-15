import pandas as pd
import numpy as np
import scipy as sp
from scipy import stats
from scipy import sparse
import pgenlib as pg
import os
import re
import warnings
import gdreg
import psutil
from sys import getsizeof
import time


################################################################################
################################## Computation #################################
################################################################################
def sample_mvn(mat_cov, random_seed=None, verbose=False):
    """
    Generate samples from zero-mean multi-variate Gaussian distribution.
    Only use positive eigenvalues from the covariance matrix.
    """
    if random_seed is not None:
        np.random.seed(random_seed)
    mat_cov = np.array(mat_cov)
    v_w, mat_v = np.linalg.eigh(mat_cov)
    ind_select = v_w > 0
    if verbose:
        print("Call: gdreg.util.sample_mvn")
        tot_var = np.absolute(v_w).sum()
        pos_var = v_w[v_w > 0].sum()
        print(
            "    dim=%d, tot_abs_eig=%0.2g, percent_pos=%0.1f%%, n_pos_eig/n_neg_eig=%d/%d"
            % (
                mat_cov.shape[0],
                tot_var,
                pos_var / tot_var * 100,
                (v_w > 0).sum(),
                (v_w < 0).sum(),
            )
        )
        print("    Top 5 eigenvalues : %s" % np.sort(v_w)[::-1][:5])
    v_w = v_w[ind_select]
    v_w = np.sqrt(v_w)
    mat_v = mat_v[:, ind_select]
    mat_u = mat_v.dot(np.diag(v_w))
    v_mvn = mat_u.dot(np.random.randn(mat_u.shape[1]))
    v_mvn = v_mvn.astype(np.float32)
    return v_mvn


def reg(mat_Y, mat_X):
    """Regress mat_Y against mat_X

    Parameters
    ----------
    mat_Y : np.ndarray
        Response variable of shape (n_sample, n_response).
    mat_X : np.ndarray
        Covariates of shape (n_sample, n_covariate).

    Returns
    -------
    mat_coef : np.ndarray
        Regression coefficients of shape (n_covariate, n_response).
    """

    if sparse.issparse(mat_X):
        mat_X = mat_X.toarray()
    else:
        mat_X = np.array(mat_X)
    if len(mat_X.shape) == 1:
        mat_X = mat_X.reshape([-1, 1])

    if sparse.issparse(mat_Y):
        mat_Y = mat_Y.toarray()
    else:
        mat_Y = np.array(mat_Y)
    if len(mat_Y.shape) == 1:
        mat_Y = mat_Y.reshape([-1, 1])

    n_sample = mat_Y.shape[0]
    mat_xtx = np.dot(mat_X.T, mat_X) / n_sample
    mat_xty = np.dot(mat_X.T, mat_Y) / n_sample
    mat_coef = np.linalg.solve(mat_xtx, mat_xty)

    if mat_coef.shape[1] == 1:
        mat_coef = mat_coef.reshape([-1])

    return mat_coef


def meta_analysis(effects, se, method="random", weights=None):
    """Random effect meta analysis"""
    # From Omer Weissbrod
    assert method in ["fixed", "random"]
    d = effects
    variances = se**2

    # compute random-effects variance tau2
    vwts = 1.0 / variances
    fixedsumm = vwts.dot(d) / vwts.sum()
    Q = np.sum(((d - fixedsumm) ** 2) / variances)
    df = len(d) - 1
    tau2 = np.maximum(0, (Q - df) / (vwts.sum() - vwts.dot(vwts) / vwts.sum()))

    # defing weights
    if weights is None:
        if method == "fixed":
            wt = 1.0 / variances
        else:
            wt = 1.0 / (variances + tau2)
    else:
        wt = weights

    # compute summtest
    summ = wt.dot(d) / wt.sum()
    if method == "fixed":
        varsum = np.sum(wt * wt * variances) / (np.sum(wt) ** 2)
    else:
        varsum = np.sum(wt * wt * (variances + tau2)) / (np.sum(wt) ** 2)
    ###summtest = summ / np.sqrt(varsum)

    summary = summ
    se_summary = np.sqrt(varsum)

    return summary, se_summary


def ldspec_meta(res_tau_list, res_rho_list, term, row, weights=None):
    """
    LD-SPEC meta-analysis. 
    - tau,h2,h2s,rho,cov,ecov : mean,se,p obtained by meta-analyzing `term / h2`
    - h2_enrich : mean,se obtained by meta-analyzing `term`, p obtained by meta-analyzing 
        `h2(c) / M(c) - [h2(maf) - h2(c)] / [M(maf) - M(c)]` (not implemented)
    - h2_shrink : mean,se obtained by meta-analyzing `term`, p obtained by meta-analyzing 
        `h2(c) - h2s(c)`
    - cor : mean,se obtained by meta-analyzing `term`, p obtained by meta-analyzing `cov`
    - ecor : mean,se obtained by meta-analyzing `term`, p obtained by meta-analyzing `ecov`
    
    Parameters
    ----------
    res_tau_list : list of DataFrame
        List of LD-SPEC tau DataFrame across traits.
    res_rho_list : list of DataFrame
        List of LD-SPEC rho DataFrame across traits.
    term : string
        One of `tau`, `h2`, `h2s`, `h2_enrich`, `h2_shrink`, `rho`, `cov`, `cor`, 
        `ecov`, `ecor`.
    row : string
        One row (annotation) in the corresponding result DataFrame.
    weights : np.ndarray (list-like), default=None
        Meta-analysis weights. If not specified, weights are equal to 1 / h2_se^2. 
    
    Returns
    -------
    meta_mean,meta_se,meta_p : floats
        Meta-analyzed mean, se, p    
        
    Todo
    ----
    - h2_enrich p not implemented
    """
    v_h2 = np.array([x.loc['AN:all', 'h2'] for x in res_tau_list], dtype=np.float32)
    v_h2_se = np.array([x.loc['AN:all', 'h2_se'] for x in res_tau_list], dtype=np.float32)
#     if weights is None: # weights
#         weights = 1 / v_h2_se**2
    
    if term in ['tau', 'h2', 'h2s']:
        # mean,se,p: meta-analyzing `term/h2` 
        v_mean = np.array([x.loc[row, term] for x in res_tau_list], dtype=np.float32)
        v_se = np.array([x.loc[row, '%s_se'%term] for x in res_tau_list], dtype=np.float32)
        meta_mean,meta_se = meta_analysis(v_mean/v_h2, v_se/v_h2, weights=weights)
        meta_p = zsc2pval(meta_mean / meta_se) 
    elif term in ['rho', 'cov', 'ecov']:
        # mean,se,p: meta-analyzing `term/h2` 
        v_mean = np.array([x.loc[row, term] for x in res_rho_list], dtype=np.float32)
        v_se = np.array([x.loc[row, '%s_se'%term] for x in res_rho_list], dtype=np.float32)
        meta_mean,meta_se = meta_analysis(v_mean/v_h2, v_se/v_h2, weights=weights)
        meta_p = zsc2pval(meta_mean / meta_se) 
    elif term in ['h2_enrich']:
        # mean,se,p: meta-analyzing `term` 
        v_mean = np.array([x.loc[row, term] for x in res_tau_list], dtype=np.float32)
        v_se = np.array([x.loc[row, '%s_se'%term] for x in res_tau_list], dtype=np.float32)
        meta_mean,meta_se = meta_analysis(v_mean, v_se, weights=weights)
        meta_p = zsc2pval( (meta_mean-1) / meta_se) 
    elif term == 'h2_shrink':
        # mean,se: meta-analyzing `term` 
        v_mean = np.array([x.loc[row, term] for x in res_tau_list], dtype=np.float32)
        v_se = np.array([x.loc[row, '%s_se'%term] for x in res_tau_list], dtype=np.float32)
        meta_mean,meta_se = meta_analysis(v_mean, v_se, weights=weights)
        # p: meta-analyzing `h2(c) - h2s(c)`
        mean_list = [x.loc[row, 'h2'] - x.loc[row, 'h2s'] for x in res_tau_list]
        z_list = [zsc2pval(x.loc[row, 'h2_shrink_p'], option='two-sided') for x in res_tau_list]
        se_list = [x/y for x,y in zip(mean_list, z_list)]
        temp_mean,temp_se = meta_analysis(np.array(mean_list)/v_h2, np.array(se_list)/v_h2, weights=weights)
        meta_p = zsc2pval(temp_mean / temp_se) 
    elif term == 'cor':
        # mean,se: meta-analyzing `term` 
        v_mean = np.array([x.loc[row, term] for x in res_rho_list], dtype=np.float32)
        v_se = np.array([x.loc[row, '%s_se'%term] for x in res_rho_list], dtype=np.float32)
        meta_mean,meta_se = meta_analysis(v_mean, v_se, weights=weights)
        # p: meta-analyzing `cov`
        v_mean = np.array([x.loc[row, 'cov'] for x in res_rho_list], dtype=np.float32)
        v_se = np.array([x.loc[row, 'cov_se'] for x in res_rho_list], dtype=np.float32)
        temp_mean,temp_se = meta_analysis(v_mean/v_h2, v_se/v_h2, weights=weights)
        meta_p = zsc2pval(temp_mean / temp_se) 
    elif term == 'ecor':
        # mean,se: meta-analyzing `term` 
        v_mean = np.array([x.loc[row, term] for x in res_rho_list], dtype=np.float32)
        v_se = np.array([x.loc[row, '%s_se'%term] for x in res_rho_list], dtype=np.float32)
        meta_mean,meta_se = meta_analysis(v_mean, v_se, weights=weights)
        # p: meta-analyzing `ecov`
        v_mean = np.array([x.loc[row, 'ecov'] for x in res_rho_list], dtype=np.float32)
        v_se = np.array([x.loc[row, 'ecov_se'] for x in res_rho_list], dtype=np.float32)
        temp_mean,temp_se = meta_analysis(v_mean/v_h2, v_se/v_h2, weights=weights)
        meta_p = zsc2pval(temp_mean / temp_se) 
    else:
        raise ValueError("term '%s' not supported" % term)
        
    return meta_mean,meta_se,meta_p


def zsc2pval(zsc, option="two-sided"):
    """
    Convert z-score to p-value. Accurate up to `zsc=36` and `pval=4.2e-284`.
    """
    #     return 1 - sp.stats.norm.cdf(zsc)
    if option == "one-sided":
        return sp.stats.norm.cdf(-zsc)  # This is more accurate
    if option == "two-sided":
        return sp.stats.norm.cdf(-np.absolute(zsc)) * 2
    
    
def pval2zsc(pval, option="two-sided"):
    """
    Convert p-value to a positive z-score. Accurate up to `zsc=36` and `pval=4.2e-284`.
    """
    if option == "one-sided":
        return -sp.stats.norm.ppf(pval)
    if option == "two-sided":
        return -sp.stats.norm.ppf(pval/2)

# def pval2zsc(pval):
#     """
#     Convert one-sided p-value to z-score. Accurate up to `zsc=36` and `pval=4.2e-284`.
#     """
#     return -sp.stats.norm.ppf(pval


# def ldspec_meta(res_list, term, row):
#     """
#     Meta-analysis for LD-SPEC estimates.
#     """
#     if term in ['tau', 'h2', 'h2_enrich', 'rho', 'cov', 'ecov']: # use the term itself as test statistics
#         mean_list = [x.loc[row, term] for x in res_list]
#         se_list = [x.loc[row, '%s_se'%term] for x in res_list]
#     elif term == 'h2_shrink': # use h2-h2s as test statistics
#         mean_list = [x.loc[row, 'h2'] - x.loc[row, 'h2s'] for x in res_list]
#         z_list = [zsc2pval(x.loc[row, 'h2_shrink_p'], option='two-sided') for x in res_list]
#         se_list = [x/y for x,y in zip(mean_list, z_list)]
#     elif term == 'cor': # use `cov` as test statistics
#         mean_list = [x.loc[row, 'cov'] for x in res_list]
#         se_list = [x.loc[row, 'cov_se'] for x in res_list]
#     elif term == 'ecor': # use `ecov` as test statistics
#         mean_list = [x.loc[row, 'ecov'] for x in res_list]
#         se_list = [x.loc[row, 'ecov_se'] for x in res_list]
#     else:
#         raise ValueError("term '%s' not supported" % term)

#     meta_mean,meta_se = meta_analysis(np.array(mean_list), np.array(se_list), weights=np.ones(len(mean_list)))
# #     meta_mean,meta_se = meta_analysis(np.array(mean_list), np.array(se_list))
#     meta_z = meta_mean / meta_se
#     meta_p = zsc2pval(meta_z) 
#     return meta_mean,meta_se,meta_z,meta_p


################################################################################
################################# LDSPEC model #################################
################################################################################
def pannot_to_csr(v_gene, sym_non_pAN="non-pAN", flag_matS=True):
    """Convert .pannot to a sparse indicator matrix

    Parameters
    ----------
    v_gene : np.ndarray
        Gene annotation.
    sym_non_pAN : str
        Symbol for SNPs not on the gene annotation.
    flag_matS : bool
        If true, return mat_S; else, return mat_G.

    Returns
    -------
    mat_S : sp.sparse.csr(dtype=bool)
        SNP-gene indicator matrix of shape (n_snp, n_gene); column corresponding
        to `sym_non_pAN` excluded.
    mat_G : sp.sparse.csr(dtype=bool)
        SNP-pair annotation indicator matrix of shape (n_snp, n_snp). Diagonal elements
        are zero.
    """
    v_gene = np.array(v_gene)
    n_snp = v_gene.shape[0]
    ind_pAN = v_gene != sym_non_pAN

    dic_gene = {}
    for i, gene in zip(np.arange(n_snp)[ind_pAN], v_gene[ind_pAN]):
        dic_gene.setdefault(gene, []).append(i)

    row_ind = []
    col_ind = []
    for i, gene in enumerate(dic_gene):
        row_ind.extend(dic_gene[gene])
        col_ind.extend([i] * len(dic_gene[gene]))

    mat_S = sp.sparse.csr_matrix(
        ([True] * len(row_ind), (row_ind, col_ind)),
        shape=[n_snp, len(dic_gene)],
        dtype=bool,
    )

    if flag_matS:
        return mat_S
    else:
        mat_G = mat_S.dot(mat_S.T)
        mat_G.setdiag(
            False
        )  # This step is taking a bit longer time. Use mat_S for HP computing
        #         mat_G = sp.sparse.csr_matrix(mat_G)
        mat_G.eliminate_zeros()
        return mat_G


def pannot_hr_to_csr(v_snp, snp_pair_list):
    """Convert .pannot_hr to a sparse indicator matrix

    Parameters
    ----------
    v_snp : np.ndarray (list-like)
        List of SNPs.
    snp_pair_list : list of lists
        List of SNP pairs.

    Returns
    -------
    mat_G : sp.sparse.csr(dtype=bool)
        SNP-pair annotation indicator matrix of shape (n_snp, n_snp). Diagonal
        elements are zero.
    """
    n_snp = len(v_snp)
    dic_snp = {snp: i for i, snp in enumerate(v_snp)}

    ind_set = set()
    for s1, s2 in snp_pair_list:
        ind_set.add((dic_snp[s1], dic_snp[s2]))
        ind_set.add((dic_snp[s2], dic_snp[s1]))

    row_ind = []
    col_ind = []
    for i, j in ind_set:
        row_ind.append(i)
        col_ind.append(j)

    mat_G = sp.sparse.csr_matrix(
        ([True] * len(row_ind), (row_ind, col_ind)), shape=[n_snp, n_snp], dtype=bool
    )
    return mat_G


def pannot_hr_to_csr_block(v_snp, snp_pair_list, block_size=int(1e6), verbose=False):
    """Convert .pannot_hr to a sparse indicator matrix by blocks to avoid large set operations

    Parameters
    ----------
    v_snp : np.ndarray (list-like)
        List of SNPs.
    snp_pair_list : list of lists
        List of SNP pairs.
    block_size : int

    Returns
    -------
    mat_G : sp.sparse.csr(dtype=bool)
        SNP-pair annotation indicator matrix of shape (n_snp, n_snp). Diagonal
        elements are zero.
    """
    n_snp = len(v_snp)
    dic_snp = {snp: i for i, snp in enumerate(v_snp)}
    mat_G = None

    n_block = np.ceil(len(snp_pair_list) / block_size).astype(int)

    for i_block in range(n_block):
        if verbose:
            print("Block %d/%d" % (i_block, n_block))
        ind_set = set()
        ind_s = i_block * block_size
        ind_e = min((i_block + 1) * block_size, len(snp_pair_list))
        for s1, s2 in snp_pair_list[ind_s:ind_e]:
            ind_set.add((dic_snp[s1], dic_snp[s2]))
            ind_set.add((dic_snp[s2], dic_snp[s1]))

        row_ind = []
        col_ind = []
        for i, j in ind_set:
            row_ind.append(i)
            col_ind.append(j)

        temp_mat_G = sp.sparse.csr_matrix(
            ([True] * len(row_ind), (row_ind, col_ind)),
            shape=[n_snp, n_snp],
            dtype=bool,
        )

        if mat_G is None:
            mat_G = temp_mat_G.copy()
        else:
            mat_G = mat_G + temp_mat_G
    return mat_G


################################################################################
################################## Read/write ##################################
################################################################################
def from_filepattern(fp, sub="@"):
    """Find all files following a given pattern

    Parameters
    ----------
    fp : str
        File pattern. `sub` can represent 'a-zA-Z0-9' of any length.
    sub : str
        Symbol used to define `fp`.

    Returns
    -------
    flist : list
        List of files following the filepattern.
    """

    sep = os.path.sep
    fpath = sep.join(fp.split(sep)[:-1])
    fname = fp.split(sep)[-1].replace(sub, "([_a-zA-Z0-9]+)")
    fname = "^" + fname + "$"

    flist = [
        "%s%s%s" % (fpath, sep, x)
        for x in os.listdir(fpath)
        if re.match(r"%s" % fname, x)
    ]
    return flist


def read_pgen(
    plink_prefix, pgen_file=None, psam_file=None, pvar_file=None, afreq_file=None
):
    """
    Load PLINK2 .pgen data

    Parameters
    ----------
    plink_prefix : str
        PLINK2 pfile prefix (for .pgen, .psam, .pvar, and .afreq files).

    pgen_file : str
        PLINK2 .pgen file.

    psam_file : str
        PLINK2 .psam file.

    pvar_file : str
        PLINK2 .pvar file.

    afreq_file : str
        PLINK2 .afreq file.

    Returns
    -------
    dic_data : dic

    TODO
    ----
    Make afreq_file optional

    """

    if pgen_file is None:
        pgen_file = plink_prefix + ".pgen"

    if psam_file is None:
        psam_file = plink_prefix + ".psam"

    if pvar_file is None:
        pvar_file = plink_prefix + ".pvar"

    if afreq_file is None:
        afreq_file = plink_prefix + ".afreq"

    assert os.path.exists(pgen_file), "pgen_file does not exist: %s" % pgen_file
    assert os.path.exists(psam_file), "psam_file does not exist: %s" % psam_file
    assert os.path.exists(pvar_file), "pvar_file does not exist: %s" % pvar_file
    assert os.path.exists(afreq_file), "afreq_file does not exist: %s" % afreq_file

    dic_data = {}

    # .psam
    n_sharp_line = get_n_sharp_line(psam_file)
    dic_data["psam"] = pd.read_csv(
        psam_file, skiprows=n_sharp_line - 1, delim_whitespace=True
    )
    dic_data["psam"].columns = update_columns(dic_data["psam"].columns)

    # .pvar
    n_sharp_line = get_n_sharp_line(pvar_file)
    dic_data["pvar"] = pd.read_csv(
        pvar_file, skiprows=n_sharp_line - 1, delim_whitespace=True
    )
    dic_data["pvar"].columns = update_columns(dic_data["pvar"].columns)

    # .pgen using pgenlib
    dic_data["pgen"] = pgen_file
    n_sample, n_snp = dic_data["psam"].shape[0], dic_data["pvar"].shape[0]
    mat_X = read_geno(pgen_file, 0, 5, n_sample, n_snp)

    # .afreq
    n_sharp_line = get_n_sharp_line(afreq_file)
    dic_data["afreq"] = pd.read_csv(
        afreq_file, skiprows=n_sharp_line - 1, delim_whitespace=True
    )
    dic_data["afreq"].columns = update_columns(dic_data["afreq"].columns)

    return dic_data


def read_geno(pgen_file, ind_start, ind_end, n_sample=None, n_snp=None):
    """Read genotype from .pgen using pg.PgenReader

    Parameters
    ----------
    pgen_file : str
        PLINK2 .pgen file

    ind_start : int
        Start index

    ind_end : int
        End index

    n_sample : int
        Number of samples in the .pgen file (obtained from .psam file).

    n_snp : int
        Total number of SNPs in the .pgen file (obtained from .pvar file).

    Returns
    -------
    mat_X : np.array(np.int8)
        Genotype matrix

    """

    with pg.PgenReader(bytes(pgen_file, encoding="utf8"), n_sample, n_snp) as reader:
        if n_sample is None:
            n_sample = reader.get_raw_sample_ct()
        mat_X = np.empty([ind_end - ind_start, n_sample], np.int8)
        reader.read_range(ind_start, ind_end, mat_X)
    return mat_X


def get_n_sharp_line(fpath):
    """Return number of header lines starting with '#'"""
    n_line = 0
    with open(fpath, "r") as f:
        for line in f:
            if line.startswith("#"):
                n_line += 1
            else:
                break
    return n_line


def update_columns(col_list):
    """
    Update columns to a set of standardized names

    Parameters
    ----------
    col_list : list
        List of columns.

    Returns
    -------
    col_list_new : list
        List of updated columns.

    """

    col_list = list(col_list)
    col_list_new = []
    for col in col_list:

        col_new = col

        if col in ["#FID", "fid", "FID", "Fid"]:
            col_new = "FID"
        if col in ["iid", "IID", "Iid"]:
            col_new = "IID"
        if col in ["chr", "CHR", "chrom", "CHROM", "#CHROM"]:
            col_new = "CHR"
        if col in ["snp", "SNP", "id", "ID"]:
            col_new = "SNP"
        if col in ["pos", "POS", "position", "POSITION", "bp", "BP"]:
            col_new = "BP"
        if col in ["cm", "CM", "pos_cm", "POS_CM"]:
            col_new = "CM"
        if col in ["gene", "GENE"]:
            col_new = "GENE"
        if col in ["n", "N", "n_sample", "N_SAMPLE"]:
            col_new = "N"
        if col in ["z", "Z", "z_score", "Z_SCORE"]:
            col_new = "Z"
        if col in ["A1", "a1", "allele1", "ALLELE1", "ALT"]:
            col_new = "ALT"
        if col in ["A2", "a2", "allele2", "ALLELE2", "REF"]:
            col_new = "REF"
        if col in ["maf", "MAF", "ALT_FREQS"]:
            col_new = "MAF"

        if col in ["EFF", "eff", "effect", "BETA", "beta"]:
            col_new = "EFFECT"

        col_list_new.append(col_new)

    return col_list_new


def write_annot(df_annot, fpath):
    """
    Write annotation files, in `.tsv.gz` format.

        - `.annot.gz` : single-SNP annotation, same as in LDSC.
        - `.pannot.gz` : SNP-pair annotation. All SNP pairs with
            the same `pAN` value belong to annotation `pAN`.
        - `.pannot_hr.gz` : high-resolution SNP-pair annotation. Each
            row is one annotated SNP pair.

    Example of `.annot.gz` file:

        CHR   SNP    BP    CM    AN:AN1    AN:AN2
        10    10:92981:C:T    92981    0    1   0.5
        10    10:93015:G:A    93015    0    5   0.1

    Example of `.pannot.gz` file (note `non-pAN` is used for SNPs outside annotation):
        CHR   SNP    BP    pAN:pAN1
        1    1:69761:A:T    69761    non-pAN
        10    10:92981:C:T    92981    Gene1
        10    10:93015:G:A    93015    Gene1

    Example for `.pannot_hr.gz` file:
        CHR   SNP    BP    pCHR   pSNP    pBP    pAN:pAN1
        10    10:92981:C:T    92981    10    10:93015:G:A    93015    1

    Parameters
    ----------
    df_annot : pd.DataFrame
        Annotations.
    fpath : str
        Output file.
    """

    # Check fpath
    err_msg = (
        "fpath should end with one of ['.annot.gz', '.pannot.gz', '.pannot_hr.gz'] : %s"
        % fpath
    )
    assert (
        fpath.endswith(".annot.gz")
        | fpath.endswith("pannot.gz")
        | fpath.endswith("pannot_hr.gz")
    ), err_msg

    # Check columns
    df_annot.columns = update_columns(df_annot.columns)
    assert "CHR" in df_annot.columns, "'CHR' missing from df_annot.columns"
    assert "SNP" in df_annot.columns, "'SNP' missing from df_annot.columns"
    assert "BP" in df_annot.columns, "'BP' missing from df_annot.columns"

    if fpath.endswith("pannot_hr.gz"):
        assert "pCHR" in df_annot.columns, "'pCHR' missing from df_annot.columns"
        assert "pSNP" in df_annot.columns, "'pSNP' missing from df_annot.columns"
        assert "pBP" in df_annot.columns, "'pBP' missing from df_annot.columns"

    # Write .annot.gz file
    if fpath.endswith(".annot.gz"):
        if "CM" not in df_annot.columns:
            print("'CM' missing from df_annot.columns, add 'CM' column with 0")
            df_annot["CM"] = 0
        col_list = ["CHR", "SNP", "BP", "CM"] + sorted(
            [x for x in df_annot.columns if x.startswith("AN:")]
        )
        df_annot[col_list].to_csv(
            fpath, index=False, header=True, sep="\t", compression="gzip"
        )

    # Write .pannot.gz file
    if fpath.endswith(".pannot.gz"):
        col_list = ["CHR", "SNP", "BP"] + sorted(
            [x for x in df_annot.columns if x.startswith("pAN:")]
        )
        df_annot[col_list].to_csv(
            fpath, index=False, header=True, sep="\t", compression="gzip"
        )

    # Write .pannot_hr.gz file
    if fpath.endswith(".pannot_hr.gz"):
        col_list = ["CHR", "SNP", "BP", "pCHR", "pSNP", "pBP"] + sorted(
            [x for x in df_annot.columns if x.startswith("pAN:")]
        )
        df_annot[col_list].to_csv(
            fpath, index=False, header=True, sep="\t", compression="gzip"
        )
    return


def read_annot(fpath, nrows=None):
    """
    Read annotation files, in `.tsv.gz` format.

        - `.annot.gz` : single-SNP annotation, same as in LDSC.
        - `.pannot.gz` : SNP-pair annotation. All SNP pairs with
            the same `pAN` value belong to annotation `pAN`.
        - `.pannot_hr.gz` : high-resolution SNP-pair annotation. Each
            row is one annotated SNP pair.

    Parameters
    ----------
    fpath : str
        Annotation file path.

    Returns
    -------
    df_annot : pd.DataFrame
        Annotation file.
    """

    # Check fpath
    err_msg = (
        "fpath should end with one of ['.annot.gz', '.pannot.gz', '.pannot_hr.gz'] : %s"
        % fpath
    )
    assert (
        fpath.endswith(".annot.gz")
        | fpath.endswith("pannot.gz")
        | fpath.endswith("pannot_hr.gz")
    ), err_msg

    # Read .annot.gz file
    if fpath.endswith(".annot.gz"):
        df_annot = pd.read_csv(fpath, sep="\t", nrows=nrows, compression="gzip")

    # Read .pannot.gz file
    if fpath.endswith(".pannot.gz"):
        df_annot = pd.read_csv(fpath, sep="\t", nrows=nrows, compression="gzip")

    # Read .pannot_hr.gz file
    if fpath.endswith(".pannot_hr.gz"):
        df_annot = pd.read_csv(fpath, sep="\t", nrows=nrows, compression="gzip")

    # Check columns
    df_annot.columns = update_columns(df_annot.columns)
    assert "CHR" in df_annot.columns, "'CHR' missing from df_annot.columns"
    assert "SNP" in df_annot.columns, "'SNP' missing from df_annot.columns"
    assert "BP" in df_annot.columns, "'BP' missing from df_annot.columns"

    if fpath.endswith("pannot_hr.gz"):
        assert "pCHR" in df_annot.columns, "'pCHR' missing from df_annot.columns"
        assert "pSNP" in df_annot.columns, "'pSNP' missing from df_annot.columns"
        assert "pBP" in df_annot.columns, "'pBP' missing from df_annot.columns"

    # Check `non-pAN` in `pAN` for .pannot.gz
    if fpath.endswith("pannot.gz"):
        for pAN in [x for x in df_annot.columns if x.startswith("pAN:")]:
            avg_non_pAN = (df_annot[pAN] == "non-pAN").mean()
            if avg_non_pAN < 0.05:
                warnings.warn(
                    "Annotation '%s' contains less than %0.1f%% non-pAN values"
                    % (pAN, avg_non_pAN * 100)
                )

    return df_annot


def get_annot_name_from_file(annot_file):
    """
    Get annotation name from file path. annot_file must end with
    `.annot.gz` or `.pannot_mat.npz`.
    """
    annot_name = annot_file.split(os.path.sep)[-1]
    # Suffix
    if annot_name.endswith(".annot.gz"):
        annot_type = "AN:"
        annot_name = annot_name.replace(".annot.gz", "")
    elif annot_name.endswith(".pannot_mat.npz"):
        annot_type = "pAN:"
        annot_name = annot_name.replace(".pannot_mat.npz", "")
    else:
        raise ValueError("annot_file must end with '.annot.gz' or '.pannot_mat.npz'")
    # CHR indicator
    annot_name = annot_name.replace("CHR@", "")
    annot_name = annot_name.replace("chr@", "")
    annot_name = annot_name.replace("C@", "")
    annot_name = annot_name.replace("c@", "")
    # Formatting
    annot_name = annot_name.replace(".", "_")
    annot_name = annot_name.replace("__", "_")
    annot_name = annot_name.strip("._,:")
    annot_name = annot_type + annot_name
    return annot_name


def write_pannot_mat(snp_pair_list, snp_list, prefix_out):
    """
    Write SNP-pair annotation matrix '.pannot_mat.npz' with respect to .pvar SNPs.

    Parameters
    ----------
    snp_pair_list : list of lists
        List of SNP pairs
    snp_list : list
        List of SNPs
    prefix_out : str
        Output prefix
    """
    mat_G = pannot_hr_to_csr_block(snp_list, snp_pair_list, verbose=False)
    sp.sparse.save_npz(prefix_out + ".pannot_mat", mat_G)


def read_pannot_mat(fpath):
    """
    Read SNP-pair annotation matrix '.pannot_mat.npz' files.

    Parameters
    ----------
    fpath : str
        Annotation file path.

    Returns
    -------
    snp_pair_list : pd.DataFrame
        List of SNP pairs.
    pAN : str
        Annotation name.
    """

    # Check fpath
    err_msg = "fpath should end with '.pannot_mat.npz'"
    assert fpath.endswith(".pannot_mat.npz"), err_msg

    mat_pannot = sp.sparse.load_npz(fpath)
    return mat_pannot


def read_ld(fpath):
    """
    Read LD files.

        - `_fullld.npy` : full_ld matrix, np.array(dtype=np.float32)
        - `_ld.npz` : ld matrix with SNPs in 10MB window, sp.sparse.csc_matrix(dtype=np.float32)

    Parameters
    ----------
    fpath: str
        LD file path.

    Returns
    -------
    mat_ld : np.array(dtype=np.float32) or sp.sparse.csc_matrix(dtype=np.float32)
        LD matrix of dimension (n_ref_snp, n_snp)
    dic_range : dict

        - dic_range['chr'] : chromosome
        - dic_range['start'] : start position
        - dic_range['end'] : end position
        - dic_range['chr_ref'] : reference chromosome list (List)
    """

    # Check fpath
    err_msg = "fpath should end with one of ['_fullld.npy', '_ld.npz'] : %s" % fpath
    assert fpath.endswith("_fullld.npy") | fpath.endswith("_ld.npz"), err_msg

    if fpath.endswith("_fullld.npy"):
        mat_ld = np.load(fpath)
        temp_str = [x for x in fpath.split(".") if x.endswith("_fullld")][0]
        dic_range = parse_snp_range(temp_str)

    if fpath.endswith("_ld.npz"):
        mat_ld = sp.sparse.load_npz(fpath)
        temp_str = [x for x in fpath.split(".") if x.endswith("_ld")][0]
        dic_range = parse_snp_range(temp_str)

    return mat_ld, dic_range


def parse_snp_range(snp_range):
    """Get range of SNPs to analyze.

    Parameters
    ----------
    snp_range: str
        Example: 'c1_s0_e2000_r1'

    Returns
    -------
    dic_range : dict

        - dic_range['chr'] : chromosome
        - dic_range['start'] : start position
        - dic_range['end'] : end position
        - dic_range['chr_ref'] : reference chromosome list (List)
    """

    dic_range = {x: None for x in ["chr", "start", "end", "chr_ref"]}

    for x in snp_range.split("_"):

        if x[0] == "c":
            dic_range["chr"] = int(x.replace("c", "").strip())

        if x[0] == "s":
            dic_range["start"] = int(x.replace("s", "").strip())

        if x[0] == "e":
            dic_range["end"] = int(x.replace("e", "").strip())

        if x[0] == "r":
            temp_str = x.replace("r", "").strip()
            if temp_str == "all":
                dic_range["chr_ref"] = list(np.arange(1, 23))
            else:
                dic_range["chr_ref"] = [int(x) for x in temp_str.split(",")]

    return dic_range


# def parse_snp_range(snp_range):
#     """Get range of SNPs to analyze. Example: 'chr=1|start=0|end=500|chr_ref=2'"""

#     dic_snp_range = {x: None for x in ["chr", "start", "end", "chr_ref"]}

#     for x in snp_range.split("|"):

#         if x.startswith("chr="):
#             dic_snp_range["chr"] = int(x.replace("chr=", "").strip())

#         if x.startswith("start="):
#             dic_snp_range["start"] = int(x.replace("start=", "").strip())

#         if x.startswith("end="):
#             dic_snp_range["end"] = int(x.replace("end=", "").strip())

#         if x.startswith("chr_ref="):
#             temp_str = x.replace("chr_ref=", "").strip()
#             if temp_str == "all":
#                 dic_snp_range["chr_ref"] = list(np.arange(1, 23))
#             else:
#                 dic_snp_range["chr_ref"] = [int(x) for x in temp_str.split(",")]

#     return dic_snp_range


def get_mbin(maf):
    if maf >= 0.05:
        return "common"
    elif maf >= 0.005:
        return "lf"
    else:
        return "rare"


################################################################################
###################################### CLI #####################################
################################################################################
def get_cli_head():
    MASTHEAD = "******************************************************************************\n"
    MASTHEAD += "* Gene-level directional effect regression (GDREG)\n"
    MASTHEAD += "* Version %s\n" % gdreg.__version__
    MASTHEAD += "* Martin Jinye Zhang\n"
    MASTHEAD += "* HSPH / Broad Institute\n"
    MASTHEAD += "* MIT License\n"
    MASTHEAD += "******************************************************************************\n"
    return MASTHEAD


def get_memory():
    """Get memory (in GB) used by the process"""
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / 1024**3
    return mem


def get_sys_info(start_time):
    """Return system info"""
    info = "sys_time=%0.1fs, sys_mem=%0.2gGB" % (time.time() - start_time, get_memory())
    return info


def sizeof(a):
    return getsizeof(a) / 1024 / 1024


def sizeof_ndarray(a):
    return (a.nbytes) / 1024 / 1024


def sizeof_sparse(a):
    return (a.data.nbytes + a.indptr.nbytes + a.indices.nbytes) / 1024 / 1024
