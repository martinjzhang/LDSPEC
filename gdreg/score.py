import numpy as np
import scipy as sp
import pandas as pd
import pgenlib as pg
import time
import gdreg
import warnings


def compute_ld(
    dic_data,
    pos_tar,
    pos_ref,
    memory=512,
    verbose=False,
):
    """
    Compute the ld matrix for two sets of SNPs.

    Parameters
    ----------
    dic_data : dict
        Genotype data reader, organized, for each CHR, as

        - dic_data[CHR]['pgen'] : .pgen file path
        - dic_data[CHR]['pvar'] : .pvar pd.DataFrame
        - dic_data[CHR]['psam'] : .psam pd.DataFrame

    pos_tar,pos_ref : list of int
        Genomic range of SNPs. CHR,ind_start,ind_end=pos_tar.
    memory : int, default=128
        Memory to use (in MB).
    verbose : bool, default=False
        If to output messages.

    Returns
    -------
    mat_ld : np.array(dtype=np.float32)
        LD matrix of shape (n_snp_ref, n_snp_tar).
    """

    start_time = time.time()

    # Check sample
    CHR_tar, ind_s_tar, ind_e_tar = pos_tar
    CHR_ref, ind_s_ref, ind_e_ref = pos_ref

    if CHR_tar != CHR_ref:
        err_msg = "Samples do not match for CHR_tar=%d and CHR_ref=%d" % (
            CHR_tar,
            CHR_ref,
        )
        assert (
            dic_data[CHR_tar]["psam"]["IID"] != dic_data[CHR_ref]["psam"]["IID"]
        ).sum() == 0, err_msg

    n_sample = dic_data[CHR_tar]["psam"].shape[0]
    n_snp_tar = ind_e_tar - ind_s_tar
    n_snp_ref = ind_e_ref - ind_s_ref

    # Open PgenReader only once
    reader_tar = pg.PgenReader(bytes(dic_data[CHR_tar]["pgen"], encoding="utf8"))
    reader_ref = pg.PgenReader(bytes(dic_data[CHR_ref]["pgen"], encoding="utf8"))

    # block_size_tar, block_size_ref, block_size_sample
    mem_ld = 2 * (n_snp_ref * n_snp_tar * 4) / 1024**2
    mem_geno_ps = (
        n_sample * (1 + 0.5 * 2) / 1024**2
    )  # Sparse matrix takes half of memory
    block_size_tot = (memory - mem_ld - 20) / mem_geno_ps
    block_size_tar = int(block_size_tot * n_snp_ref / (n_snp_tar + n_snp_ref))
    block_size_ref = int(block_size_tot * n_snp_tar / (n_snp_tar + n_snp_ref))
    block_size_sample = 16383

    err_msg = "block_size_tot=%d too small, allocate at least %dMB memory" % (
        block_size_tot,
        mem_ld + mem_geno_ps * 2048 + 20,
    )
    assert block_size_tot > 2048, err_msg

    n_block_tar = np.ceil(n_snp_tar / block_size_tar).astype(int)
    n_block_ref = np.ceil(n_snp_ref / block_size_ref).astype(int)
    block_size_tar = np.ceil(n_snp_tar / n_block_tar).astype(int)
    block_size_ref = np.ceil(n_snp_ref / n_block_ref).astype(int)
    n_block_sample = np.ceil(n_sample / block_size_sample).astype(int)

    if verbose:
        print("# Call: gdreg.score.compute_ld")
        print(
            "    n_snp_tar=%d (CHR%d), n_snp_ref=%d (CHR%d), n_sample=%d"
            % (n_snp_tar, CHR_tar, n_snp_ref, CHR_ref, n_sample)
        )
        print("    memory=%dMB" % (memory))
        print("    block_size_tar=%d, n_block_tar=%d" % (block_size_tar, n_block_tar))
        print("    block_size_ref=%d, n_block_ref=%d" % (block_size_ref, n_block_ref))
        print(
            "    block_size_sample=%d, n_block_sample=%d"
            % (block_size_sample, n_block_sample)
        )

    # Compute LD matrix
    mat_ld = np.zeros([n_snp_ref, n_snp_tar], dtype=int)
    v_maf_tar, v_maf_ref = [], []

    for i_tar in range(n_block_tar):
        start_tar = ind_s_tar + i_tar * block_size_tar
        end_tar = min(ind_s_tar + (i_tar + 1) * block_size_tar, ind_e_tar)
        mat_X_tar = np.empty([end_tar - start_tar, n_sample], np.int8)
        reader_tar.read_range(start_tar, end_tar, mat_X_tar)
        mat_X_tar[mat_X_tar == -9] = 0
        v_maf_tar.extend(mat_X_tar.mean(axis=1) * 0.5)

        mat_list_tar = []
        for i_sample in range(n_block_sample):
            ind_s = i_sample * block_size_sample
            ind_e = (i_sample + 1) * block_size_sample
            mat_list_tar.append(
                sp.sparse.csr_matrix(mat_X_tar[:, ind_s:ind_e].T, dtype=np.int16)
            )

        for i_ref in range(n_block_ref):
            start_ref = ind_s_ref + i_ref * block_size_ref
            end_ref = min(ind_s_ref + (i_ref + 1) * block_size_ref, ind_e_ref)
            mat_X_ref = np.empty([end_ref - start_ref, n_sample], np.int8)
            reader_ref.read_range(start_ref, end_ref, mat_X_ref)
            mat_X_ref[mat_X_ref == -9] = 0
            if i_tar == 0:
                v_maf_ref.extend(mat_X_ref.mean(axis=1) * 0.5)

            for i_sample in range(n_block_sample):
                ind_s = i_sample * block_size_sample
                ind_e = (i_sample + 1) * block_size_sample
                temp_mat = sp.sparse.csr_matrix(
                    mat_X_ref[:, ind_s:ind_e].T, dtype=np.int16
                )
                mat_ld[
                    start_ref - ind_s_ref : end_ref - ind_s_ref,
                    start_tar - ind_s_tar : end_tar - ind_s_tar,
                ] += temp_mat.T.dot(mat_list_tar[i_sample]).toarray()

    reader_tar.close()
    reader_ref.close()

    # Normalization
    mat_ld = mat_ld.astype(np.float32) / n_sample
    v_maf_tar = np.array(v_maf_tar, dtype=np.float32)
    v_maf_ref = np.array(v_maf_ref, dtype=np.float32)

    mat_ld = mat_ld - np.outer(2 * v_maf_ref, 2 * v_maf_tar)
    temp_v1 = 1 / np.sqrt(2 * v_maf_ref * (1 - v_maf_ref))
    temp_v2 = 1 / np.sqrt(2 * v_maf_tar * (1 - v_maf_tar))
    mat_ld = mat_ld * np.outer(temp_v1, temp_v2)

    if verbose:
        print("    Completed, time=%0.1fs" % (time.time() - start_time))
    return mat_ld


def compute_score(
    dic_data,
    dic_ld,
    df_annot,
    pannot_list=[],
    pannot_hr_list=[],
    sym_non_pAN="non-pAN",
    win_size=int(1e7),
    snp_range=None,
    memory=512,
    verbose=False,
    verbose_prefix="",
):

    """
    Compute DLD scores for all annots and all pannot-annot pairs.

    Parameters
    ----------
    dic_data : dict
        Genotype data reader, organized, for each CHR, as

        - dic_data[CHR]['pgen'] : .pgen file path
        - dic_data[CHR]['pvar'] : .pvar pd.DataFrame
        - dic_data[CHR]['psam'] : .psam pd.DataFrame

    dic_ld : dict of sp.sparse.csc_matrix
        dic_mat_ld[CHR] for LD matrix of chromosome CHR. The j-th column contains
        all LDs of SNP j within a 1e7 window (5e6 on each side).
    df_annot : pd.DataFrame
        Single-SNP annotation.
    df_pannot_list : list of pd.DataFrame, default=[]
        Each element corresponds to SNP-pair annotation. Must contain
        ['CHR', 'SNP', 'BP', 'pAN:pAN1'] columns.
    df_pannot_hr_list : list of pd.DataFrame, default=[]
        Each element corresponds to a high-res SNP-pair annotation. Must contain
        ['CHR', 'SNP', 'BP', 'pCHR', 'pSNP', 'pBP', 'pAN:pAN1'] columns.
    sym_non_pAN : str, default='non-pAN'
        Symbol for SNPs not in the SNP-pair annotation.
    win_size : int, defualt=1e7
        Window size for computing LD and DLD scores.
    snp_range : list of int
        Genomic range of SNPs for computing the scores. (CHR,ind_start,ind_end).
        If provide, only only scores for this range.
    memory : int, default=128
        Memory to use (in MB).
    verbose : bool, default=False
        If to output messages.

    Returns
    -------
    df_score : pd.DataFrame
        GDREG LD and DLD scores, with columns ['CHR', 'SNP', 'BP', 'LD:AN1', 'LD:AN2', 'LD:E', 'DLD:PAN:AN1', 'DLD:PAN:AN2'].

    TODO
    ----
    - Remove 'n_sample' for computing LD scores

    """

    start_time = time.time()

    # SNP info : check consistency between dic_data and dic_ld
    CHR_list = sorted(set(dic_data))
    for CHR in CHR_list:
        if CHR not in dic_ld:
            continue
        assert dic_ld[CHR].shape[0] == dic_data[CHR]["pvar"].shape[0], (
            "dic_ld does not match dimension of dic_data for CHR%d" % CHR
        )
    n_sample = dic_data[CHR_list[0]]["psam"].shape[0]
    v_snp = []
    for CHR in CHR_list:
        v_snp.extend(dic_data[CHR]["pvar"]["SNP"])

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

    # block_size
    # TODO : determine size
    block_size = 500

    if verbose:
        print("# Call: gdreg.score.compute_score")
        temp_str = ", ".join(
            ["CHR%d (%d SNPs)" % (x, dic_data[x]['pvar'].shape[0]) for x in CHR_list]
        )
        print("    %d SNPs from %d CHRs: %s" % (len(v_snp), len(CHR_list), temp_str))
        if snp_range is not None:
            print("    Range: chr=%d, start=%d, end=%d " % snp_range)
        print("    Single-SNP annots : %s" % ", ".join(AN_list))
        print(
            "    SNP-pair annots : %s"
            % ", ".join(pAN_list + ["%s (hr)" % x for x in pAN_hr_list])
        )
        print("    win_size=%0.1fMB, memory=%dMB" % (win_size / 1e6, memory))

    # Compute score
    df_score = None
    for CHR in CHR_list:
        
        if snp_range is not None:
            if CHR != snp_range[0]:
                continue

        dic_score = {"CHR": [], "SNP": [], "BP": [], "E": []}
        dic_score.update({"LD:%s" % x: [] for x in AN_list})
        dic_score.update({"DLD:%s" % x: [] for x in pAN_list})
        dic_score.update({"DLD:%s" % x: [] for x in pAN_hr_list})

        v_snp_chr = dic_data[CHR]["pvar"]["SNP"].values
        v_bp_chr = dic_data[CHR]["pvar"]["BP"].values
        n_snp_chr = v_snp_chr.shape[0]
        n_block_chr = np.ceil(n_snp_chr / block_size).astype(int)

        for i_block in range(n_block_chr):
            ind_s = i_block * block_size
            ind_e = min((i_block + 1) * block_size, n_snp_chr)
            ind_s_ref = np.searchsorted(
                v_bp_chr, v_bp_chr[ind_s] - win_size / 2, side="left"
            )
            ind_e_ref = np.searchsorted(
                v_bp_chr, v_bp_chr[ind_e - 1] + win_size / 2, side="right"
            )

            if snp_range is not None:
                if (ind_s > snp_range[2]) | (ind_e < snp_range[1]):
                    continue

            # Basic info
            dic_score["CHR"].extend(dic_data[CHR]["pvar"]["CHR"].values[ind_s:ind_e])
            dic_score["SNP"].extend(dic_data[CHR]["pvar"]["SNP"].values[ind_s:ind_e])
            dic_score["BP"].extend(dic_data[CHR]["pvar"]["BP"].values[ind_s:ind_e])

            mat_ld_block = dic_ld[CHR][:, ind_s:ind_e].toarray()[
                ind_s_ref:ind_e_ref, :
            ]  # (n_ref, n_tar)
            v_snp_ref_block = list(v_snp_chr[ind_s_ref:ind_e_ref])

            # E score : r_ii
            n_dif = ind_s - ind_s_ref
            dic_score["E"].extend(
                [mat_ld_block[x + n_dif, x] for x in range(ind_e - ind_s)]
            )

            # LD score
            for AN in AN_list:
                v_annot = [
                    dic_annot[AN][x] if x in dic_annot[AN] else 0
                    for x in v_snp_ref_block
                ]
                v_annot = np.array(v_annot, dtype=np.float32)
                v_annot_outside = [
                    dic_annot[AN][x] if x in dic_annot[AN] else 0
                    for x in set(v_snp) - set(v_snp_ref_block)
                ]
                v_score = (mat_ld_block.T**2 * v_annot).sum(axis=1)
                v_score += np.sum(v_annot_outside) / n_sample
                dic_score["LD:%s" % AN].extend(v_score)

            # DLD score
            # New DLD score for non-interacting pannot
            for pAN in pAN_list:
                v_pAN = [
                    dic_pannot[pAN][x] if x in dic_pannot[pAN] else sym_non_pAN
                    for x in v_snp_ref_block
                ]
                mat_S = gdreg.util.pannot_to_csr(v_pAN)
                mat_G = mat_S.dot(mat_S.T).toarray()
                np.fill_diagonal(mat_G, 0)
                v_score = (mat_G.dot(mat_ld_block) * mat_ld_block).sum(axis=0)
                dic_score["DLD:%s" % pAN].extend(v_score)

            for pAN in pAN_hr_list:
                snp_set = set(v_snp_ref_block)
                snp_pair_list = [
                    x
                    for x in dic_pannot_hr[pAN]
                    if (x[0] in snp_set) & (x[1] in snp_set)
                ]
                mat_G = gdreg.util.pannot_hr_to_csr(v_snp_ref_block, snp_pair_list)
                v_score = (mat_G.dot(mat_ld_block) * mat_ld_block).sum(axis=0)
                dic_score["DLD:%s" % pAN].extend(v_score)

        temp_df = pd.DataFrame(dic_score)
        if df_score is None:
            df_score = temp_df.copy()
        else:
            df_score = pd.concat([df_score, temp_df], axis=0)

    if snp_range is not None:
        CHR, START, END = snp_range
        snp_list = dic_data[CHR]["pvar"]["SNP"].values[START:END]
        df_score = df_score.loc[df_score["SNP"].isin(snp_list)]

    if verbose:
        print(
            verbose_prefix + "    Completed, time=%0.1fs" % (time.time() - start_time)
        )

    return df_score


# def compute_dld_score(mat_ld, mat_G, v_annot, v_ps_sd):
#     """
#     Compute directional LD score:

#         `s_i = r_i^T (mat_G * 0.5(v_annot 1^T + 1^T v_annot ) * (v_ps_sd v_ps_sd^T)) r_i`

#     where `*` means hadamard product; dot product is omitted. `r_i` means LD between
#     SNP i and all other reference SNPs.


#     Parameters
#     ----------
#     mat_ld : np.ndarray
#         LD matrix of shape (n_snp_ref, n_snp)
#     mat_G : np.ndarray
#         SNP-pair indicator of shape (n_snp_ref, n_snp_ref)
#     v_annot : np.ndarray
#         Annotation of shape (n_snp_ref,)
#     v_ps_sd : np.ndarray
#         Per-SNP causal effect SD of shape (n_snp_ref,)

#     Returns
#     -------
#     v_score : np.array(dtype=np.float32)
#         DLD score of shape (n_snp,)
#     """
#     if sp.sparse.issparse(mat_G):
#         mat_sigma = mat_G.multiply(v_ps_sd).T.multiply(v_annot * v_ps_sd).toarray()
#     else:
#         mat_sigma = (mat_G * v_ps_sd).T * (v_annot * v_ps_sd)
#     v_score = (mat_sigma.dot(mat_ld) * mat_ld).sum(axis=0)
#     return v_score
