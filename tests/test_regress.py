import pytest
import numpy as np
import scipy as sp
import pandas as pd
from scipy import sparse
import ldspec
import os


def test_reg_bjn():
    """
    Test ldspec.regress.reg_bjn
    """

    n_sample = 100
    n_rep = 20

    v_beta = np.array([1, -1])
    mat_coef = np.zeros([n_rep, 2])
    ten_cov = np.zeros([n_rep, 2, 2])
    for i_rep in range(n_rep):
        np.random.seed(i_rep)
        mat_X = np.random.randn(n_sample, 2)
        v_y = mat_X.dot(v_beta) + np.random.randn(n_sample) * 1
        dic_block = ldspec.regress.get_block(
            pd.DataFrame(data={"CHR": [1] * mat_X.shape[0]}), n_block=10
        )
        dic_res = ldspec.regress.reg_bjn(v_y, mat_X, dic_block, verbose=False)
        mat_coef[i_rep] = dic_res["coef"]
        ten_cov[i_rep] = dic_res["coef_cov"]

    esti_mean = mat_coef.mean(axis=0)
    rel_abs_err = np.absolute(v_beta - esti_mean).sum() / np.absolute(v_beta).sum()
    err_msg = "mean: true=%s, esti=%s, rel_abs_err=%0.3e" % (
        str(v_beta),
        str(esti_mean),
        rel_abs_err,
    )
    print(err_msg)
    assert rel_abs_err < 0.05, err_msg

    true_se = np.sqrt(np.diag(np.cov(mat_coef.T)))
    esti_se = np.sqrt(np.diag(ten_cov.mean(axis=0)))
    mean_se_ratio = np.mean(esti_se / true_se)
    err_msg = "se: true=%s, esti=%s\nse_ratio=%s, mean_se_ratio=%0.3e" % (
        str(true_se),
        str(esti_se),
        str(esti_se / true_se),
        mean_se_ratio,
    )
    print(err_msg)
    assert mean_se_ratio > 0.9, err_msg

    return
