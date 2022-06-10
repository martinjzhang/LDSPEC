import pandas as pd
import numpy as np
import scipy as sp
from scipy import sparse
import pgenlib as pg
import os
import warnings
import gdreg
import psutil
from sys import getsizeof
import time


################################################################################
################################## Computation #################################
################################################################################
def sample_mvn(mat_cov, random_seed=0, verbose=False):
    """
    Generate samples from zero-mean multi-variate Gaussian distribution.
    Only use positive eigenvalues from the covariance matrix.
    """
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
    v_mvn[np.absolute(v_mvn) < 1e-8] = 0
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


################################################################################
################################## GDReg model #################################
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
    df_annot: pd.DataFrame
        Annotations.

    fpath: str
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


def read_annot(fpath):
    """
    Read annotation files, in `.tsv.gz` format.

        - `.annot.gz` : single-SNP annotation, same as in LDSC.
        - `.pannot.gz` : SNP-pair annotation. All SNP pairs with
            the same `pAN` value belong to annotation `pAN`.
        - `.pannot_hr.gz` : high-resolution SNP-pair annotation. Each
            row is one annotated SNP pair.

    Parameters
    ----------
    fpath: str
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
        df_annot = pd.read_csv(fpath, sep="\t", compression="gzip")

    # Read .pannot.gz file
    if fpath.endswith(".pannot.gz"):
        df_annot = pd.read_csv(fpath, sep="\t", compression="gzip")

    # Read .pannot_hr.gz file
    if fpath.endswith(".pannot_hr.gz"):
        df_annot = pd.read_csv(fpath, sep="\t", compression="gzip")

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


def parse_snp_range(snp_range):
    """Get range of SNPs to analyze. Example: 'chr=1|start=0|end=500|chr_ref=2'"""

    dic_snp_range = {x: None for x in ["chr", "start", "end", "chr_ref"]}

    for x in snp_range.split("|"):

        if x.startswith("chr="):
            dic_snp_range["chr"] = int(x.replace("chr=", "").strip())

        if x.startswith("start="):
            dic_snp_range["start"] = int(x.replace("start=", "").strip())

        if x.startswith("end="):
            dic_snp_range["end"] = int(x.replace("end=", "").strip())

        if x.startswith("chr_ref="):
            temp_str = x.replace("chr_ref=", "").strip()
            if temp_str == "all":
                dic_snp_range["chr_ref"] = list(np.arange(1, 23))
            else:
                dic_snp_range["chr_ref"] = [int(x) for x in temp_str.split(",")]

    return dic_snp_range


################################################################################
###################################### CLI #####################################
################################################################################
def get_cli_head():
    MASTHEAD = "******************************************************************************\n"
    MASTHEAD += "* Gene-level directional effect regression (GDReg)\n"
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
