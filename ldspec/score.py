import numpy as np
import scipy as sp
import pandas as pd
import pgenlib as pg
import time
import ldspec
import warnings


def compute_ld(
    dic_data,
    pos_tar,
    pos_ref,
    verbose=False,
):
    """
    Compute the ld matrix for two sets of SNPs.

    Parameters
    ----------
    dic_data : dict
        Genotype data reader. Must contain all CHRs.
        - dic_data[CHR]['pgen'] : .pgen file path
        - dic_data[CHR]['pvar'] : .pvar pd.DataFrame
        - dic_data[CHR]['psam'] : .psam pd.DataFrame
    pos_tar,pos_ref : list of int
        Genomic range of SNPs of format [CHR, ind_start, ind_end].

    Returns
    -------
    mat_ld : np.array(dtype=np.float32)
        LD matrix of shape (n_snp_ref, n_snp_tar).
    """

    start_time = time.time()

    CHR_tar, ind_s_tar, ind_e_tar = pos_tar
    CHR_ref, ind_s_ref, ind_e_ref = pos_ref

    # Check sample
    if CHR_tar != CHR_ref:
        err_msg = "Samples don't match for CHR_tar=%d and CHR_ref=%d" % (
            CHR_tar,
            CHR_ref,
        )
        assert (
            dic_data[CHR_tar]["psam"]["IID"] != dic_data[CHR_ref]["psam"]["IID"]
        ).sum() == 0, err_msg

    n_sample = dic_data[CHR_tar]["psam"].shape[0]
    n_snp_tar = ind_e_tar - ind_s_tar
    n_snp_ref = ind_e_ref - ind_s_ref

    # block_size_tar, block_size_ref, block_size_sample
    block_size_tar = n_snp_tar
    block_size_ref = 1000
    block_size_sample = 10000  # 32767/4

    n_block_tar = np.ceil(n_snp_tar / block_size_tar).astype(int)
    n_block_ref = np.ceil(n_snp_ref / block_size_ref).astype(int)
    block_size_tar = np.ceil(n_snp_tar / n_block_tar).astype(int)
    block_size_ref = np.ceil(n_snp_ref / n_block_ref).astype(int)
    n_block_sample = np.ceil(n_sample / block_size_sample).astype(int)

    if verbose:
        print("# Call: ldspec.score.compute_ld")
        print(
            "    n_snp_tar=%d (CHR%d), n_snp_ref=%d (CHR%d), n_sample=%d"
            % (n_snp_tar, CHR_tar, n_snp_ref, CHR_ref, n_sample)
        )
        print("    block_size_tar=%d, n_block_tar=%d" % (block_size_tar, n_block_tar))
        print("    block_size_ref=%d, n_block_ref=%d" % (block_size_ref, n_block_ref))
        print(
            "    block_size_sample=%d, n_block_sample=%d"
            % (block_size_sample, n_block_sample)
        )

    # Open PgenReader only once
    reader_tar = pg.PgenReader(bytes(dic_data[CHR_tar]["pgen"], encoding="utf8"))
    reader_ref = pg.PgenReader(bytes(dic_data[CHR_ref]["pgen"], encoding="utf8"))

    # Compute LD matrix
    mat_ld = np.zeros([n_snp_ref, n_snp_tar], dtype=float)
    v_maf_tar, v_maf_ref = [], []

    for i_tar in range(n_block_tar):
        start_tar = ind_s_tar + i_tar * block_size_tar
        end_tar = min(ind_s_tar + (i_tar + 1) * block_size_tar, ind_e_tar)
        mat_X_tar = np.empty([end_tar - start_tar, n_sample], np.int8)
        reader_tar.read_range(start_tar, end_tar, mat_X_tar)
        v_missing = (mat_X_tar == -9).sum(axis=1)  # Imputation by mean genotype
        v_maf_tar_block = (
            0.5 * (mat_X_tar.sum(axis=1) + 9 * v_missing) / (n_sample - v_missing)
        )
        v_maf_tar.extend(v_maf_tar_block)  # np.float64

        mat_list_tar = []
        for i_sample in range(n_block_sample):
            ind_s = i_sample * block_size_sample
            ind_e = (i_sample + 1) * block_size_sample
            temp_mat = sp.sparse.csr_matrix(
                mat_X_tar[:, ind_s:ind_e].T, dtype=np.float32
            )
            temp_ind = temp_mat == -9
            temp_mat += 9 * temp_ind + temp_ind.multiply(v_maf_tar_block * 2)
            mat_list_tar.append(temp_mat)

        for i_ref in range(n_block_ref):
            start_ref = ind_s_ref + i_ref * block_size_ref
            end_ref = min(ind_s_ref + (i_ref + 1) * block_size_ref, ind_e_ref)
            mat_X_ref = np.empty([end_ref - start_ref, n_sample], np.int8)
            reader_ref.read_range(start_ref, end_ref, mat_X_ref)
            v_missing = (mat_X_ref == -9).sum(axis=1)  # Imputation by mean genotype
            v_maf_ref_block = (
                0.5 * (mat_X_ref.sum(axis=1) + 9 * v_missing) / (n_sample - v_missing)
            )
            if i_tar == 0:
                v_maf_ref.extend(v_maf_ref_block)  # np.float64

            for i_sample in range(n_block_sample):
                ind_s = i_sample * block_size_sample
                ind_e = (i_sample + 1) * block_size_sample
                temp_mat = sp.sparse.csr_matrix(
                    mat_X_ref[:, ind_s:ind_e].T, dtype=np.float32
                )
                temp_ind = temp_mat == -9
                temp_mat += 9 * temp_ind + temp_ind.multiply(v_maf_ref_block * 2)
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
    mat_ld[np.isnan(mat_ld)] = 0
    mat_ld[np.isinf(mat_ld)] = 0
    # make diagonal values 1, assuming n_snp_tar<n_snp_ref
    if (ind_s_tar >= ind_s_ref) & (ind_e_tar <= ind_e_ref):
        for i_tar in range(n_snp_tar):
            i_ref = i_tar + ind_s_tar - ind_s_ref
            if mat_ld[i_ref, i_tar] == 0:
                mat_ld[i_ref, i_tar] = 1

    if verbose:
        print("    Completed, time=%0.1fs" % (time.time() - start_time))
    return mat_ld


def compute_avgr(
    dic_pannot_path,
    dic_ld_path,
    verbose=False,
):
    """
    Compute average LD for SNP-pair annotations.

    Parameters
    ----------
    dic_pannot_path : dict of dict of strs
        File path for SNP-pair annotations. `dic_pannot_path[annot_name][CHR]` contains the
        `.pannot_mat.npz` file path for annotation `annot_name` and `CHR`. Dimension of the
        sparse matrix must match `dic_data[CHR]['pvar']`.
    dic_ld_path : dict of dict of strs
        File path for LD matrices of shape (n_snp_ref, n_snp_tar), in csc format.
        dic_ld_path[CHR] contains the list of LD files for `CHR`. Dimension of the
        sparse matrix should match `dic_data[CHR]['pvar']`.

    Returns
    -------
    dic_avgr : dic, default={}
        dic_avgr[pAN] contains the average LD across all pairs in pAN.
    """
    start_time = time.time()
    np.random.seed(0)
    pAN_list = list(dic_pannot_path)
    dic_sumr = {x: 0 for x in pAN_list}
    dic_n = {x: 0 for x in pAN_list}
    for CHR in dic_ld_path:
        if verbose:
            print("CHR%d (%d LD files)" % (CHR, len(dic_ld_path[CHR])))
        dic_mat_G_chr = {}
        for pAN in pAN_list:
            dic_mat_G_chr[pAN] = ldspec.util.read_pannot_mat(dic_pannot_path[pAN][CHR])

        p_thresold = 5 / len(dic_ld_path[CHR])
        for i, ld_file in enumerate(dic_ld_path[CHR]):
            if np.random.rand(1)[0] > p_thresold:
                continue
            mat_ld, dic_range = ldspec.util.read_ld(ld_file)
            mat_ld.data[np.isnan(mat_ld.data)] = 0
            for pAN in pAN_list:
                temp_mat_G = dic_mat_G_chr[pAN][
                    dic_range["start"] : dic_range["end"], :
                ].T
                dic_sumr[pAN] += temp_mat_G.multiply(mat_ld).sum()
                dic_n[pAN] += temp_mat_G.sum()
            if verbose:
                print(
                    "    LD file %d/%d, time=%0.1fs"
                    % (i, len(dic_ld_path[CHR]), time.time() - start_time)
                )
    dic_avgr = {x: dic_sumr[x] / dic_n[x] for x in pAN_list}
    if verbose:
        print("    Completed, time=%0.1fs" % (time.time() - start_time))
    return dic_avgr


def compute_score(
    dic_data,
    dic_ld,
    dic_annot_path={},
    dic_pannot_path={},
    snp_range=None,
    flag_cross_term=False,
    win_size=int(1e6),
    verbose=False,
):

    """
    Compute LD/DLD scores for annots/pannots. If `snp_range` is given, compute scores for
    SNPs in the range. Otherwise, compute scores for all SNPs on CHRs which with both
    `dic_data[CHR]` and `dic_ld[CHR]`.

    Parameters
    ----------
    dic_data : dict
        Genotype data reader. Must contain all CHRs.
        - dic_data[CHR]['pgen'] : .pgen file path
        - dic_data[CHR]['pvar'] : .pvar pd.DataFrame
        - dic_data[CHR]['psam'] : .psam pd.DataFrame
    dic_ld : dict of sp.sparse.csc_matrix(dtype=np.float32)
        dic_ld[CHR] contains the signed LD matrix (r) for CHR with shape (n_snp_chr, n_snp_chr).
        The j-th column of dic_ld[CHR] contains the signed LD between SNP j and all other SNPs on the CHR.
        Only values within `win_size` of SNP j will be used (`win_size/2` on each side).
        If `snp_range` is given, columns specified by `snp_range` will be used and other columns can
        be filled with zero.
    dic_annot_path : dict of dict of strs
        File path for single-SNP annotations. `dic_annot_path[annot_name][CHR]` contains the
        `.annot.gz` file path for `annot_name` and `CHR`.
    dic_pannot_path : dict of dict of strs
        File path for SNP-pair annotations. `dic_pannot_path[annot_name][CHR]` contains the
        `.pannot_mat.npz` file path for annotation `annot_name` and `CHR`. Dimension of the
        sparse matrix must match `dic_data[CHR]['pvar']`.
    snp_range : list of int
        SNP range of the format `(CHR, ind_start, ind_end)`.
        If provided, only compute scores for SNPs in this range.
    flag_cross_term : bool, default=False
        If True, also compute scores for cross terms (Z_i Z_j) for SNP pairs i,j within 10000 SNPs
        and covered by at least one pannot. Not tested.
    win_size : int, defualt=1e6
        Window size.

    Returns
    -------
    df_score : pd.DataFrame
        LD and DLD scores, with columns ['CHR', 'SNP', 'BP', 'E', 'LD:AN1', 'LD:AN2', ...,
        'DLD:pAN1', 'DLD:pAN2', ...].

    TODO
    ----

    """

    start_time = time.time()

    # SNP info
    CHR_list = sorted(set(dic_data) & set(dic_ld))
    for CHR in CHR_list:
        # Check dic_ld
        assert (
            dic_ld[CHR].shape[0] == dic_ld[CHR].shape[1]
        ), "dic_ld[%d] is not a squared matrix: %s" % (CHR, str(dic_ld[CHR].shape))
        assert (
            dic_ld[CHR].shape[0] == dic_data[CHR]["pvar"].shape[0]
        ), "dic_ld[%d] does not match dimension of dic_data[%d]" % (CHR, CHR)
        # Check dic_annot_path
        for annot_name in dic_annot_path:
            assert (
                CHR in dic_annot_path[annot_name]
            ), "`dic_annot_path['%s'][%d] missing`" % (annot_name, CHR)
        # Check dic_pannot_path
        for annot_name in dic_pannot_path:
            assert (
                CHR in dic_pannot_path[annot_name]
            ), "`dic_pannot_path['%s'][%d] missing`" % (annot_name, CHR)

    # Annotation
    AN_list = []
    for annot_name in dic_annot_path:
        temp_df = ldspec.util.read_annot(
            dic_annot_path[annot_name][CHR_list[0]], nrows=5
        )
        AN_list.extend([x for x in temp_df if x.startswith("AN:")])
    pAN_list = list(dic_pannot_path)

    if verbose:
        print("# Call: ldspec.score.compute_score")
        temp_str = ", ".join(
            ["CHR%d (%d SNPs)" % (x, dic_data[x]["pvar"].shape[0]) for x in CHR_list]
        )
        print("    CHRs in `dic_data` and `dic_ld`: %s" % (temp_str))
        if snp_range is not None:
            print("    Range: chr=%d, start=%d, end=%d " % snp_range)
        print("    Annots: %s" % ", ".join(AN_list))
        print("    Pannots: %s" % ", ".join(pAN_list))
        print("    win_size=%0.1fMB" % (win_size / 1e6))

    # Compute score
    block_size = 1000
    df_score_list = []
    for CHR in CHR_list:  # Iterate over all CHRs, process one block a time
        if snp_range is not None:
            if CHR != snp_range[0]:
                continue

        dic_score = {"CHR": [], "SNP": [], "BP": [], "E": []}
        dic_score.update({"LD:%s" % x: [] for x in AN_list})
        dic_score.update({"DLD:%s" % x: [] for x in pAN_list})

        v_snp_chr = dic_data[CHR]["pvar"]["SNP"].values
        v_bp_chr = dic_data[CHR]["pvar"]["BP"].values
        n_snp_chr = v_snp_chr.shape[0]
        n_block_chr = np.ceil(n_snp_chr / block_size).astype(int)

        # Read annots and pannots for CHR
        dic_annot_chr = {}
        for annot_name in dic_annot_path:
            temp_df = ldspec.util.read_annot(dic_annot_path[annot_name][CHR])
            for AN in [x for x in temp_df if x.startswith("AN:")]:
                dic_annot_chr[AN] = {
                    x: y for x, y in zip(temp_df["SNP"], temp_df[AN]) if y != 0
                }
        dic_mat_G_chr = {}
        for pAN in pAN_list:
            dic_mat_G_chr[pAN] = ldspec.util.read_pannot_mat(dic_pannot_path[pAN][CHR])

        for i_block in range(n_block_chr):
            ind_s = i_block * block_size
            ind_e = min((i_block + 1) * block_size, n_snp_chr)
            if snp_range is not None:
                if (ind_s >= snp_range[2]) | (ind_e <= snp_range[1]):
                    continue
            # Possible bug: dic_ld doesn't contain all values
            # between ind_s_ref and ind_e_ref due to block mismatch.
            # This should be minor though.
            ind_s_ref = np.searchsorted(
                v_bp_chr, v_bp_chr[ind_s] - win_size / 2, side="left"
            )
            ind_e_ref = np.searchsorted(
                v_bp_chr, v_bp_chr[ind_e - 1] + win_size / 2, side="right"
            )

            # Basic info
            dic_score["CHR"].extend([CHR] * (ind_e - ind_s))
            dic_score["SNP"].extend(v_snp_chr[ind_s:ind_e])
            dic_score["BP"].extend(v_bp_chr[ind_s:ind_e])

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
                    dic_annot_chr[AN][x] if x in dic_annot_chr[AN] else 0
                    for x in v_snp_ref_block
                ]
                v_annot = np.array(v_annot, dtype=np.float32)
                v_score = (mat_ld_block.T**2 * v_annot).sum(axis=1)
                dic_score["LD:%s" % AN].extend(v_score)

            # DLD score
            for pAN in pAN_list:
                temp_mat_G = dic_mat_G_chr[pAN][
                    ind_s_ref:ind_e_ref, ind_s_ref:ind_e_ref
                ].copy()
                v_score = (temp_mat_G.dot(mat_ld_block) * mat_ld_block).sum(axis=0)
                dic_score["DLD:%s" % pAN].extend(v_score)

            # Cross terms (not tested)
            if flag_cross_term:
                # ind_pair
                v_snp_tar_block = list(v_snp_chr[ind_s:ind_e])
                v_bp_tar_block = list(v_bp_chr[ind_s:ind_e])
                n_snp_tar_block = len(v_snp_tar_block)
                mat_pair = np.zeros([n_snp_tar_block, n_snp_tar_block], dtype=bool)
                for pAN in pAN_list:
                    mat_pair = (
                        mat_pair
                        | dic_mat_G_chr[pAN][ind_s:ind_e, ind_s:ind_e].toarray()
                    )

                ind_select = mat_pair[np.triu_indices(n_snp_tar_block, k=1)]
                temp_ = np.triu_indices(n_snp_tar_block, k=1)
                ind_pair = [
                    (x, y) for x, y in zip(temp_[0][ind_select], temp_[1][ind_select])
                ]
                n_pair = len(ind_pair)

                # Basic info
                dic_score["CHR"].extend([CHR] * n_pair)
                dic_score["SNP"].extend(
                    [
                        "%s|%s" % (v_snp_tar_block[x], v_snp_tar_block[y])
                        for x, y in ind_pair
                    ]
                )
                dic_score["BP"].extend(
                    [
                        int(0.5 * (v_bp_tar_block[x] + v_bp_tar_block[y]))
                        for x, y in ind_pair
                    ]
                )

                # E score : r_ij
                n_dif = ind_s - ind_s_ref
                row_list = [n_dif + x[0] for x in ind_pair]
                col_list = [x[1] for x in ind_pair]
                dic_score["E"].extend(mat_ld_block[(row_list, col_list)])

                # LD score : \sum_k r_ik r_jk a_ck + r_ij * n_pair_outside / n_sample
                # The second term not computed for computational efficiency
                for AN in AN_list:
                    v_annot = [
                        dic_annot_chr[AN][x] if x in dic_annot_chr[AN] else 0
                        for x in v_snp_ref_block
                    ]
                    v_annot = np.array(v_annot, dtype=np.float32)
                    mat_score = (mat_ld_block.T * v_annot).dot(mat_ld_block)
                    row_list = [x[0] for x in ind_pair]
                    col_list = [x[1] for x in ind_pair]
                    dic_score["LD:%s" % AN].extend(mat_score[(row_list, col_list)])

                # DLD score : \sum_k,k' r_ik r_jk' G_kk'
                for pAN in pAN_list:
                    temp_mat_G = dic_mat_G_chr[pAN][
                        ind_s_ref:ind_e_ref, ind_s_ref:ind_e_ref
                    ].copy()
                    mat_score = temp_mat_G.dot(mat_ld_block).T.dot(mat_ld_block)
                    row_list = [x[0] for x in ind_pair]
                    col_list = [x[1] for x in ind_pair]
                    dic_score["DLD:%s" % pAN].extend(mat_score[(row_list, col_list)])

        df_score_list.append(pd.DataFrame(dic_score))
    df_score = pd.concat(df_score_list, axis=0)

    if snp_range is not None:  # Retain only SNPs in `snp_range`
        CHR, START, END = snp_range
        snp_set = set(dic_data[CHR]["pvar"]["SNP"].values[START:END])
        ind_select = [
            (x.split("|")[0] in snp_set) | (x.split("|")[-1] in snp_set)
            for x in df_score["SNP"]
        ]
        df_score = df_score.loc[ind_select]

    if verbose:
        print("    Completed, time=%0.1fs" % (time.time() - start_time))

    return df_score
