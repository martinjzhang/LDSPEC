import pandas as pd
import numpy as np
import scipy as sp
import os
from os.path import join
import re
import time
import argparse
import pickle

# in-house tools
import gdreg


"""
Job desecription
----------------

get_snp_block : generate SNP blocks based on sample size and memory
    - Input : --job | --pgen_file | --prefix_out
    - Output : list of snp ranges (in the format of "snp_range")

compute_ld : compute LD matrix.
    - Input : --job | --pgen_file | --prefix_out | --snp_range | [--memory] ｜ [--random_seed]
        | [--flag_full_ld]
    - Output : LD matrix between a set of SNPs and all other SNPs on the same chromosome.
    
compute_score : compute LD and DLD scores.
    - Input : --job | --pgen_file | --ld_file | --annot_file | --prefix_out | [--memory]
    ｜ [--random_seed]
    - Output : LD and DLD scores.
    
regress : infer parameters \tau and \rho.
    - Input : --job | --pgen_file | --score_file | --sumstats_file | --annot_file | --prefix_out 
        | [--memory]
    - Output : GDREG result.
    
TODO
----

"""


def main(args):
    sys_start_time = time.time()

    ###########################################################################################
    ######                                    Parse Options                              ######
    ###########################################################################################

    JOB = args.job
    PGEN_FILE = args.pgen_file
    LD_FILE = args.ld_file
    SCORE_FILE = args.score_file
    SUMSTATS_FILE = args.sumstats_file
    ANNOT_FILE = args.annot_file
    PREFIX_OUT = args.prefix_out
    MEMORY = args.memory
    RANDOM_SEED = args.random_seed
    SNP_RANGE = args.snp_range
    FLAG_FULL_LD = args.flag_full_ld

    # Parse and check arguments
    LEGAL_JOB_LIST = ["get_snp_block", "compute_ld", "compute_score", "regress"]
    err_msg = "# run_gdreg: --job=%s not supported" % JOB
    assert JOB in LEGAL_JOB_LIST, err_msg

    if JOB in ["compute_score"]:
        assert LD_FILE is not None, "--ld_file required for --job=%s" % JOB
    if JOB in ["regress"]:
        assert SCORE_FILE is not None, "--score_file required for --job=%s" % JOB
    if JOB in ["regress"]:
        assert SUMSTATS_FILE is not None, "--sumstats_file required for --job=%s" % JOB
    if JOB in ["compute_score", "regress"]:
        assert ANNOT_FILE is not None, "--annot_file required for --job=%s" % JOB
    if JOB in ["compute_ld"]:
        assert SNP_RANGE is not None, "--snp_range required for --job=%s" % JOB
        DIC_RANGE = gdreg.util.parse_snp_range(SNP_RANGE)

    # Print input options
    header = gdreg.util.get_cli_head()
    header += "Call: run_gdreg.py \\\n"
    header += "--job %s\\\n" % JOB
    header += "--pgen_file %s\\\n" % PGEN_FILE
    header += "--ld_file %s\\\n" % LD_FILE
    header += "--score_file %s\\\n" % SCORE_FILE
    header += "--sumstats_file %s\\\n" % SUMSTATS_FILE
    header += "--annot_file %s\\\n" % ANNOT_FILE
    header += "--prefix_out %s\\\n" % PREFIX_OUT
    header += "--snp_range %s\\\n" % SNP_RANGE
    header += "--memory %d\\\n" % MEMORY
    header += "--random_seed %d\\\n" % RANDOM_SEED
    header += "--flag_full_ld %s\n" % FLAG_FULL_LD
    print(header)

    ###########################################################################################
    ######                                   Data Loading                                ######
    ###########################################################################################
    # Load --pgen_file
    if JOB in ["get_snp_block", "compute_ld", "compute_score", "regress"]:
        print("# Loading --pgen_file")
        dic_data = {}
        if "@" not in PGEN_FILE:
            temp_dic = gdreg.util.read_pgen(PGEN_FILE)
            dic_data[temp_dic["pvar"]["CHR"][0]] = temp_dic.copy()
        else:
            for CHR in range(1, 23):
                if os.path.exists(PGEN_FILE.replace("@", "%s" % CHR) + ".pgen"):
                    dic_data[CHR] = gdreg.util.read_pgen(
                        PGEN_FILE.replace("@", "%s" % CHR)
                    )

        for CHR in dic_data:
            n_sample = dic_data[CHR]["psam"].shape[0]
            n_snp = dic_data[CHR]["pvar"].shape[0]
            mat_X = gdreg.util.read_geno(
                dic_data[CHR]["pgen"], 0, 50, n_sample=None, n_snp=None
            )
            sparsity = (mat_X != 0).mean()
            print(
                "    CHR%2d: %d samples, %d SNPs, %0.1f%% non-zeros for first 50 SNPs"
                % (CHR, n_sample, n_snp, sparsity * 100)
            )
        print("    " + gdreg.util.get_sys_info(sys_start_time))

    # Load --ld_file
    if JOB in ["compute_score"]:
        print("# Loading --ld_file")
        assert os.path.exists(LD_FILE), "--ld_file does not exist"
        mat_ld, dic_range = gdreg.util.read_ld(LD_FILE)
        if dic_range["chr_ref"] is None:
            dic_range["chr_ref"] = dic_range["chr"]
        err_msg = "n_snp=%d, mismatch with --pgen_file" % mat_ld.shape[0]
        assert mat_ld.shape[0] == dic_data[dic_range["chr"]]["pvar"].shape[0], err_msg
        print(
            "    chr=%d, start=%d, end=%d, chr_ref=%d"
            % (
                dic_range["chr"],
                dic_range["start"],
                dic_range["end"],
                dic_range["chr_ref"],
            )
        )
        print("    n_snp=%d, n_snp_ref=%d" % (mat_ld.shape[1], mat_ld.shape[0]))
        print("    LD info loaded, matching --pgen_file")
        print("    " + gdreg.util.get_sys_info(sys_start_time))

    # Load --score_file
    if JOB in ["regress"]:
        print("# Loading --score_file")
        flist = sorted(gdreg.util.from_filepattern(SCORE_FILE))
        print("    find %d score files" % len(flist))
        df_score = None
        for fpath in flist:
            temp_df = pd.read_csv(fpath, sep="\t", index_col=None)

            if df_score is None:
                df_score = temp_df.copy()
            else:
                df_score = pd.concat([df_score, temp_df], axis=0)

        df_score.sort_values(["CHR", "BP"], inplace=True)
        LD_list = [x for x in df_score if x.startswith("LD:")]
        DLD_list = [x for x in df_score if x.startswith("DLD:")]

        print(
            "    score file loaded for %d SNPs, %d LD scores, %d DLD scores"
            % (df_score.shape[0], len(LD_list), len(DLD_list))
        )
        print("    " + gdreg.util.get_sys_info(sys_start_time))

    # Load --sumstats_file
    if JOB in ["regress"]:
        print("# Loading --sumstats_file")
        df_sumstats = pd.read_csv(SUMSTATS_FILE, sep="\t", index_col=None)
        print("    .sumstats.gz loaded, %d SNPs" % df_sumstats.shape[0])
        print("    " + gdreg.util.get_sys_info(sys_start_time))

    # Load --annot_file
    if JOB in ["compute_score", "regress"]:
        print("# Loading --annot_file")
        df_annot = None
        pannot_list = []
        pannot_hr_list = []
        for annot_file in ANNOT_FILE.split(","):
            err_msg = "--annot_file missing : '%s'" % annot_file
            assert os.path.exists(annot_file), err_msg
            temp_df = gdreg.util.read_annot(annot_file)

            if annot_file.endswith(".annot.gz"):
                temp_df.index = temp_df["SNP"]
                if df_annot is None:
                    df_annot = temp_df.copy()
                else:
                    col_list = [x for x in temp_df if x.startswith("AN:")]
                    df_annot = df_annot.join(temp_df[col_list])
            if annot_file.endswith(".pannot.gz"):
                pannot_list.append(temp_df.copy())
            if annot_file.endswith(".pannot_hr.gz"):
                pannot_hr_list.append(temp_df.copy())
        AN_list = [x for x in df_annot if x.startswith("AN:")]
        print(
            "    .annot.gz (%d SNPs and %d annots): %s"
            % (df_annot.shape[0], len(AN_list), ",".join(AN_list))
        )
        temp_list = ["%s (%d SNPs)" % (x.columns[-1], x.shape[0]) for x in pannot_list]
        print(
            "    .pannot.gz (%d pannots): %s" % (len(pannot_list), ",".join(temp_list)),
        )
        temp_list = [
            "%s (%d pairs)" % (x.columns[-1], x.shape[0]) for x in pannot_hr_list
        ]
        print(
            "    .pannot_hr.gz (%d pannots): %s"
            % (len(pannot_hr_list), ",".join(temp_list)),
        )
        print("    " + gdreg.util.get_sys_info(sys_start_time))

    ###########################################################################################
    ######                                  Computation                                  ######
    ###########################################################################################
    if JOB == "get_snp_block":
        print("# Running --job get_snp_block")
        block_size = 3000
        fout = open(PREFIX_OUT + ".snp_range.txt", "w")
        for CHR in dic_data:
            n_snp = dic_data[CHR]["pvar"].shape[0]
            n_block = np.ceil(n_snp / block_size).astype(int)
            for i in range(n_block):
                START = i * block_size
                END = min((i + 1) * block_size, n_snp)
                fout.write("c%d_s%d_e%d\n" % (CHR, START, END))
        fout.close()
        print("    " + gdreg.util.get_sys_info(sys_start_time))

    if JOB == "compute_ld":
        print("# Running --job compute_ld")
        if FLAG_FULL_LD:
            CHR, CHR_REF = DIC_RANGE["chr"], DIC_RANGE["chr_ref"][0]
            pos_tar = [CHR, 0, dic_data[CHR]["pvar"].shape[0]]
            pos_ref = [CHR_REF, 0, dic_data[CHR_REF]["pvar"].shape[0]]
            mat_ld = gdreg.score.compute_ld(
                dic_data, pos_tar, pos_ref, verbose=True, memory=MEMORY
            )
            np.save(PREFIX_OUT + ".c%s_r%s_fullld" % (CHR, CHR_REF), mat_ld)
        else:
            CHR, START, END = DIC_RANGE["chr"], DIC_RANGE["start"], DIC_RANGE["end"]
            n_snp_tar, n_snp = (END - START), dic_data[CHR]["pvar"].shape[0]
            v_bp = dic_data[CHR]["pvar"]["BP"].values

            block_size = 1000
            n_block = np.ceil(n_snp_tar / block_size).astype(int)
            mat_ld_list = []
            for i_block in range(n_block):
                print(
                    "block %d/%d %s"
                    % (i_block, n_block, gdreg.util.get_sys_info(sys_start_time))
                )
                ind_s = START + i_block * block_size
                ind_e = min(START + (i_block + 1) * block_size, END)
                ind_s_ref = np.searchsorted(v_bp, v_bp[ind_s] - 5.05e6, side="left")
                ind_s_ref = max(0, ind_s_ref - 3)
                ind_e_ref = np.searchsorted(
                    v_bp, v_bp[ind_e - 1] + 5.05e6, side="right"
                )
                ind_e_ref = min(n_snp, ind_e_ref + 3)

                pos_tar = [CHR, ind_s, ind_e]
                pos_ref = [CHR, ind_s_ref, ind_e_ref]
                mat_ld = gdreg.score.compute_ld(
                    dic_data,
                    pos_tar,
                    pos_ref,
                    verbose=True,
                    memory=MEMORY,
                )
                temp_mat = np.zeros([n_snp, ind_e - ind_s], dtype=np.float32)
                temp_mat[ind_s_ref:ind_e_ref, :] = mat_ld
                mat_ld_list.append(sp.sparse.csc_matrix(temp_mat))

            mat_ld = sp.sparse.hstack(mat_ld_list, format="csc")
            sp.sparse.save_npz(PREFIX_OUT + ".%s_ld" % SNP_RANGE, mat_ld)
        print("    " + gdreg.util.get_sys_info(sys_start_time))

    if JOB == "compute_score":
        CHR, START, END = dic_range["chr"], dic_range["start"], dic_range["end"]
        n_snp = dic_data[CHR]["pvar"].shape[0]

        # Zero pad for mat_ld
        mat_ld_list = []
        if START > 0:
            mat_ld_list.append(
                sp.sparse.csc_matrix(
                    ([0], ([0], [0])), shape=[n_snp, START], dtype=np.float32
                )
            )
        mat_ld_list.append(mat_ld)
        if END < n_snp:
            mat_ld_list.append(
                sp.sparse.csc_matrix(
                    ([0], ([0], [0])), shape=[n_snp, n_snp - END], dtype=np.float32
                )
            )
        dic_ld = {CHR: sp.sparse.hstack(mat_ld_list, format="csc")}
        snp_range = (dic_range["chr"], dic_range["start"], dic_range["end"])

        df_score = gdreg.score.compute_score(
            dic_data,
            dic_ld,
            df_annot,
            pannot_list,
            pannot_hr_list,
            verbose=True,
            win_size=1e7,
            snp_range=snp_range,
        )

        df_score.to_csv(
            PREFIX_OUT + ".c%d_s%d_e%d_score.tsv.gz" % (CHR, START, END),
            sep="\t",
            index=False,
            compression="gzip",
        )

    if JOB == "regress":
        print("# Running --job regress")

        dic_res = gdreg.regress.estimate(
            dic_data,
            df_score,
            df_sumstats,
            df_annot,
            pannot_list=pannot_list,
            pannot_hr_list=pannot_hr_list,
            n_jn_block=100,
            sym_non_pAN="non-pAN",
            verbose=True,
        )

        # Store the entire file and a summary df
        dbfile = open(PREFIX_OUT + ".pickle", "wb")
        pickle.dump(dic_res, dbfile)
        dbfile.close()
        dic_res[0]["summary"]["tau"].to_csv(
            PREFIX_OUT + ".tau.tsv", sep="\t", index=False
        )
        dic_res[1]["summary"]["tau"].to_csv(
            PREFIX_OUT + ".joint_tau.tsv", sep="\t", index=False
        )
        dic_res[1]["summary"]["rho"].to_csv(
            PREFIX_OUT + ".joint_rho.tsv", sep="\t", index=False
        )

        print("    " + gdreg.util.get_sys_info(sys_start_time))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="gdreg")

    parser.add_argument("--job", type=str, required=True, help="One of [compute_ld]")
    parser.add_argument("--pgen_file", type=str, required=True)
    parser.add_argument("--ld_file", type=str, required=False, default=None)
    parser.add_argument("--score_file", type=str, required=False, default=None)
    parser.add_argument("--sumstats_file", type=str, required=False, default=None)
    parser.add_argument("--annot_file", type=str, required=False, default=None)
    parser.add_argument("--prefix_out", type=str, required=True)
    parser.add_argument(
        "--snp_range",
        type=str,
        default=None,
        help="SNP range, e.g., `chr=1|start=0|end=500|chr_ref=2`. `chr_ref=all` means all 22 chromosomes.",
    )
    parser.add_argument("--memory", type=int, default=512)
    parser.add_argument("--random_seed", type=int, default=0)
    parser.add_argument("--flag_full_ld", type=bool, default=False)

    args = parser.parse_args()
    main(args)
