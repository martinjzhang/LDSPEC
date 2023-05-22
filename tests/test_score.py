import pytest
import numpy as np
import scipy as sp
import pandas as pd
from scipy import sparse
import ldspec
import os


def test_compute_ld():
    """
    Test ldspec.score.compute_ld using the public 1000G test data
    """

    DATA_PATH = os.path.join(ldspec.__path__[0], "data")
    PGEN_FILE = os.path.join(DATA_PATH, "data.geno.c@")
    dic_data = {
        CHR: ldspec.util.read_pgen(PGEN_FILE.replace("@", "%d" % CHR)) for CHR in [1, 2]
    }

    for pos_tar, pos_ref in [
        [[1, 2, 5], [1, 6889, 6907]],
        [[1, 1059, 1068], [1, 697, 715]],
        [[1, 85, 111], [2, 999, 1021]],
        [[2, 5998, 5999], [2, 2, 5]],
    ]:
        mat_X = ldspec.util.read_geno(
            dic_data[pos_tar[0]]["pgen"], pos_tar[1], pos_tar[2]
        )
        mat_X = mat_X.T.astype(np.float32)  # n_sample, n_snp
        for i in range(
            mat_X.shape[1]
        ):  # Mean imputation for each column, not used because no missing data
            temp_v = mat_X[:, i]
            mat_X[temp_v == -9, i] = np.mean(temp_v[temp_v != -9])
        v_maf = mat_X.mean(axis=0) * 0.5
        mat_X = (mat_X - 2 * v_maf) / np.sqrt(2 * v_maf * (1 - v_maf))

        mat_Y = ldspec.util.read_geno(
            dic_data[pos_ref[0]]["pgen"], pos_ref[1], pos_ref[2]
        )
        mat_Y = mat_Y.T.astype(np.float32)
        for i in range(mat_Y.shape[1]):
            temp_v = mat_Y[:, i]
            mat_Y[temp_v == -9, i] = np.mean(temp_v[temp_v != -9])
        v_maf = mat_Y.mean(axis=0) * 0.5
        mat_Y = (mat_Y - 2 * v_maf) / np.sqrt(2 * v_maf * (1 - v_maf))

        mat_ld_true = mat_Y.T.dot(mat_X) / mat_X.shape[0]
        mat_ld = ldspec.score.compute_ld(dic_data, pos_tar, pos_ref, verbose=False)

        abs_err = np.absolute(mat_ld_true - mat_ld).mean()
        err_msg = "pos_tar=%s, pos_ref=%s, abs_err=%0.3e" % (
            str(pos_tar),
            str(pos_ref),
            abs_err,
        )
        assert abs_err < 1e-4, err_msg
    return


def test_compute_ld_ukb():
    """
    Test ldspec.score.compute_ld using the private UKB test data containing missing genotypes
    """

    DATA_PATH = "/n/groups/price/martin/LDSPEC_data/test_data_ukb"
    if os.path.exists(DATA_PATH) is False:
        return

    PGEN_FILE = os.path.join(DATA_PATH, "data.geno.c@")
    dic_data = {
        CHR: ldspec.util.read_pgen(PGEN_FILE.replace("@", "%d" % CHR)) for CHR in [1, 2]
    }

    for pos_tar, pos_ref in [
        [[1, 15, 20], [1, 71, 80]],
        [[1, 1059, 1068], [1, 697, 715]],
        [[1, 85, 111], [2, 999, 1021]],
        [[2, 5998, 5999], [2, 2, 5]],
    ]:
        mat_X = ldspec.util.read_geno(
            dic_data[pos_tar[0]]["pgen"], pos_tar[1], pos_tar[2]
        )
        print("n_missing_X = %d" % (mat_X == -9).sum())
        mat_X = mat_X.T.astype(np.float32)  # n_sample, n_snp
        for i in range(
            mat_X.shape[1]
        ):  # Mean imputation for each column, not used because no missing data
            temp_v = mat_X[:, i]
            mat_X[temp_v == -9, i] = np.mean(temp_v[temp_v != -9])
        v_maf = mat_X.mean(axis=0) * 0.5
        mat_X = (mat_X - 2 * v_maf) / np.sqrt(2 * v_maf * (1 - v_maf))

        mat_Y = ldspec.util.read_geno(
            dic_data[pos_ref[0]]["pgen"], pos_ref[1], pos_ref[2]
        )
        print("n_missing_Y = %d" % (mat_Y == -9).sum())
        mat_Y = mat_Y.T.astype(np.float32)
        for i in range(mat_Y.shape[1]):
            temp_v = mat_Y[:, i]
            mat_Y[temp_v == -9, i] = np.mean(temp_v[temp_v != -9])
        v_maf = mat_Y.mean(axis=0) * 0.5
        mat_Y = (mat_Y - 2 * v_maf) / np.sqrt(2 * v_maf * (1 - v_maf))

        mat_ld_true = mat_Y.T.dot(mat_X) / mat_X.shape[0]
        mat_ld = ldspec.score.compute_ld(dic_data, pos_tar, pos_ref, verbose=False)

        abs_err = np.absolute(mat_ld_true - mat_ld).mean()
        err_msg = "pos_tar=%s, pos_ref=%s, abs_err=%0.3e" % (
            str(pos_tar),
            str(pos_ref),
            abs_err,
        )
        assert abs_err < 1e-4, err_msg
    return


def test_compute_avgr():
    """
    Test ldspec.score.compute_avgr
    """
    return


def test_compute_score():
    """
    Test ldspec.score.compute_score
    """
    DATA_PATH = os.path.join(ldspec.__path__[0], "data")
    PGEN_FILE = os.path.join(DATA_PATH, "data.geno.c@")
    LD_FILE = os.path.join(DATA_PATH, "ldspec_res/ldspec_s2.c2_s10000_e15822_ld.npz")
    ANNOT_FILE_LIST = [
        os.path.join(DATA_PATH, "data.c@.annot.gz"),
        os.path.join(
            DATA_PATH, "proxy_0_1000.ld_n100_p100.maf_all.c@.pannot_mat.npz"
        ),
    ]

    # Load test data
    dic_data = {
        CHR: ldspec.util.read_pgen(PGEN_FILE.replace("@", "%d" % CHR)) for CHR in [1, 2]
    }
    mat_ld, dic_range = ldspec.util.read_ld(LD_FILE)
    dic_annot_path, dic_pannot_path = {}, {}
    for annot_file in ANNOT_FILE_LIST:
        annot_name = ldspec.util.get_annot_name_from_file(annot_file)
        if annot_file.endswith(".annot.gz"):
            dic_annot_path[annot_name] = {}
            for CHR in dic_data:
                dic_annot_path[annot_name][CHR] = annot_file.replace("@", "%d" % CHR)
        if annot_file.endswith(".pannot_mat.npz"):
            dic_pannot_path[annot_name] = {}
            for CHR in dic_data:
                dic_pannot_path[annot_name][CHR] = annot_file.replace("@", "%d" % CHR)

    # Compute score based on implementation
    CHR = dic_range["chr"]
    n_snp = dic_data[CHR]["pvar"].shape[0]

    mat_ld_list = []
    if dic_range["start"] > 0:
        mat_ld_list.append(
            sp.sparse.csc_matrix(
                ([0], ([0], [0])), shape=[n_snp, dic_range["start"]], dtype=np.float32
            )
        )
    mat_ld_list.append(mat_ld)
    if dic_range["end"] < n_snp:
        mat_ld_list.append(
            sp.sparse.csc_matrix(
                ([0], ([0], [0])),
                shape=[n_snp, n_snp - dic_range["end"]],
                dtype=np.float32,
            )
        )
    dic_ld = {CHR: sp.sparse.hstack(mat_ld_list, format="csc")}
    snp_range = (dic_range["chr"], dic_range["start"], dic_range["start"] + 1499)

    df_score = ldspec.score.compute_score(
        dic_data,
        dic_ld,
        dic_annot_path=dic_annot_path,
        dic_pannot_path=dic_pannot_path,
        flag_cross_term=False,
        verbose=True,
        win_size=int(1.1e6),
        snp_range=snp_range,
    )

    # Test SNPs
    test_snp_relind = [1, 2, 10, 199, 698, 1005]  # index relative to dic_range['start']
    test_snp_absind = [
        x + dic_range["start"] for x in test_snp_relind
    ]  # absolute index
    test_snp_list = dic_data[CHR]["pvar"]["SNP"][test_snp_absind].tolist()

    # Test E score
    v_score_true = [mat_ld[x + dic_range["start"], x] for x in test_snp_relind]
    temp_dic = {x: y for x, y in zip(df_score["SNP"], df_score["E"])}
    v_score = [temp_dic[x] for x in test_snp_list]

    rel_abs_err = (
        np.absolute(np.array(v_score) - np.array(v_score_true)).sum()
        / np.absolute(v_score_true).sum()
    )
    err_msg = "E score: rel_abs_err=%0.3e\n" % (rel_abs_err)
    for i_snp, snp in enumerate(test_snp_list):
        err_msg += "%-15s score=%0.5f    true_score=%0.5f\n" % (
            snp,
            v_score[i_snp],
            v_score_true[i_snp],
        )
    print(err_msg)
    assert rel_abs_err < 1e-4, err_msg

    # Test LD score
    df_annot = ldspec.util.read_annot(dic_annot_path["AN:data"][CHR])
    for AN in ["AN:DHS_Trynka_common", "AN:all_common"]:
        v_annot = df_annot[AN].values
        v_score_true = []
        for i_snp in test_snp_relind:
            temp_v = mat_ld[:, i_snp].toarray().flatten()
            v_score_true.append((temp_v**2 * v_annot).sum())

        temp_dic = {x: y for x, y in zip(df_score["SNP"], df_score["LD:%s" % AN])}
        v_score = [temp_dic[x] for x in test_snp_list]

        rel_abs_err = (
            np.absolute(np.array(v_score) - np.array(v_score_true)).sum()
            / np.absolute(v_score_true).sum()
        )
        err_msg = "LD score for %s: rel_abs_err=%0.3e\n" % (AN, rel_abs_err)
        for i_snp, snp in enumerate(test_snp_list):
            err_msg += "%-15s score=%0.5f    true_score=%0.5f\n" % (
                snp,
                v_score[i_snp],
                v_score_true[i_snp],
            )
        print(err_msg)
        assert rel_abs_err < 1e-4, err_msg

    # Test DLD score
    pAN = "pAN:proxy_0_1000_ld_n100_p100_maf_all"
    mat_G = ldspec.util.read_pannot_mat(dic_pannot_path[pAN][CHR])
    v_score_true = []
    for i_snp in test_snp_relind:
        temp_v = mat_ld[:, i_snp].toarray().flatten()
        v_score_true.append(mat_G.dot(temp_v).dot(temp_v))

    temp_dic = {x: y for x, y in zip(df_score["SNP"], df_score["DLD:%s" % pAN])}
    v_score = [temp_dic[x] for x in test_snp_list]

    rel_abs_err = (
        np.absolute(np.array(v_score) - np.array(v_score_true)).sum()
        / np.absolute(v_score_true).sum()
    )
    err_msg = "DLD score for %s: rel_abs_err=%0.3e\n" % (pAN, rel_abs_err)
    for i_snp, snp in enumerate(test_snp_list):
        err_msg += "%-15s score=%0.5f    true_score=%0.5f\n" % (
            snp,
            v_score[i_snp],
            v_score_true[i_snp],
        )
    print(err_msg)
    assert rel_abs_err < 1e-4, err_msg
    return
