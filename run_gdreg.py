import pandas as pd
import numpy as np
import scipy as sp
import os
from os.path import join
import re
import time
import argparse
import pickle
import gdreg


"""
Job description
----------------

get_snp_block : generate SNP blocks (10,000 SNPs per line)
    - Input : --job | --pgen_file | --prefix_out
    - Output : list of snp ranges (in the format of "snp_range")

compute_ld : compute LD matrix.
    - Input : --job | --pgen_file | --prefix_out | --snp_range | [--random_seed] | [--flag_full_ld]
    - Output : LD matrix between a set of SNPs and all other SNPs on the same chromosome.
    
compute_score : compute LD and DLD scores.
    - Input : --job | --pgen_file | --ld_file | --annot_file | --prefix_out | [--flag_cross_term]
    - Output : LD and DLD scores.
    
compute_avgr : compute average LD (avgr) for each pannot. "--ld_file" should contain all LD files
    - Input : --job | --pgen_file | --prefix_out | --annot_file | --ld_file
    - Output : Avg. LD for each pannot.
    
regress : infer parameters \tau and \rho.
    - Input : --job | --pgen_file | --score_file | --sumstats_file | --annot_file | --avgr_file
    | --prefix_out | [--flag_cross_term]
    - Output : GDREG result.
    
TODO
----
- Double check file consistency 
- 
"""


def main(args):
    sys_start_time = time.time()

    ###########################################################################################
    ######                                    Parse Options                              ######
    ###########################################################################################

    JOB = args.job
    PGEN_FILE = args.pgen_file
    LD_FILE = args.ld_file
    ANNOT_FILE = args.annot_file
    SCORE_FILE = args.score_file
    SUMSTATS_FILE = args.sumstats_file
    AVGR_FILE = args.avgr_file
    PREFIX_OUT = args.prefix_out
    SNP_RANGE = args.snp_range
    FLAG_FULL_LD = args.flag_full_ld
    FLAG_CROSS_TERM = args.flag_cross_term
    FLAG_NOFIL_SNP = args.flag_nofil_snp

    # Parse and check arguments
    LEGAL_JOB_LIST = [
        "get_snp_block",
        "compute_ld",
        "compute_score",
        "compute_avgr",
        "regress",
    ]
    err_msg = "# run_gdreg: --job=%s not supported" % JOB
    assert JOB in LEGAL_JOB_LIST, err_msg

    if JOB in ["get_snp_block", "compute_ld", "compute_score", "regress"]:
        assert PGEN_FILE is not None, "--pgen_file required for --job=%s" % JOB
    if JOB in ["compute_score", "compute_avgr"]:
        assert LD_FILE is not None, "--ld_file required for --job=%s" % JOB
    if JOB in ["regress"]:
        assert SCORE_FILE is not None, "--score_file required for --job=%s" % JOB
    if JOB in ["regress"]:
        assert SUMSTATS_FILE is not None, "--sumstats_file required for --job=%s" % JOB
    if JOB in ["compute_score", "compute_avgr", "regress"]:
        assert ANNOT_FILE is not None, (
            "--annot_path_file required for --job=%s" % JOB
        )
    if JOB in ["compute_ld"]:
        assert SNP_RANGE is not None, "--snp_range required for --job=%s" % JOB
        DIC_RANGE = gdreg.util.parse_snp_range(SNP_RANGE)

    # Print input options
    header = gdreg.util.get_cli_head()
    header += "Call: run_gdreg.py \\\n"
    header += "--job %s\\\n" % JOB
    header += "--pgen_file %s\\\n" % PGEN_FILE
    header += "--ld_file %s\\\n" % LD_FILE
    header += "--annot_file %s\\\n" % ANNOT_FILE
    header += "--score_file %s\\\n" % SCORE_FILE
    header += "--sumstats_file %s\\\n" % SUMSTATS_FILE
    header += "--avgr_file %s\\\n" % AVGR_FILE
    header += "--prefix_out %s\\\n" % PREFIX_OUT
    header += "--snp_range %s\\\n" % SNP_RANGE
    header += "--flag_full_ld %s\\\n" % FLAG_FULL_LD
    header += "--flag_cross_term %s\\\n" % FLAG_CROSS_TERM
    header += "--flag_nofil_snp %s\n" % FLAG_NOFIL_SNP
    print(header)

    ###########################################################################################
    ######                                  Data Loading                                 ######
    ###########################################################################################
    # Load --pgen_file
    if JOB in [
        "get_snp_block",
        "compute_ld",
        "compute_score",
        "compute_avgr",
        "regress",
    ]:
        print("# Loading --pgen_file")
        dic_data = {}
        if "@" not in PGEN_FILE:  # Load single CHR
            temp_dic = gdreg.util.read_pgen(PGEN_FILE)
            dic_data[temp_dic["pvar"]["CHR"][0]] = temp_dic.copy()
        else:
            for CHR in range(1, 23):  # Check all 23 CHRs
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

    # Load --ld_file (lazy)
    if JOB in ["compute_avgr"]:
        print("# Loading --ld_file (dic_ld_path)")
        dic_ld_path = {x: [] for x in dic_data}
        for fpath in gdreg.util.from_filepattern(LD_FILE, sub="@"):
            temp_str = [x for x in fpath.split(".") if x.endswith("_ld")][0]
            dic_range = gdreg.util.parse_snp_range(temp_str)
            if dic_range["chr"] in dic_ld_path:  # Load only files for CHRs in dic_data
                dic_ld_path[dic_range["chr"]].append(fpath)
            else:
                print("    Skip: %s" % fpath)
        for CHR in dic_ld_path:
            print("    CHR%2d: detected %d LD files" % (CHR, len(dic_ld_path[CHR])))

    # Load --ld_file
    if JOB in ["compute_score"]:
        print("# Loading --ld_file (mat_ld)")
        assert os.path.exists(LD_FILE), "--ld_file does not exist"
        mat_ld, dic_range = gdreg.util.read_ld(LD_FILE)
        mat_ld.data[np.isnan(mat_ld.data)] = 0
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

    # Load --sumstats_file
    if JOB in ["regress"]:
        print("# Loading --sumstats_file")
        df_sumstats = pd.read_csv(SUMSTATS_FILE, sep="\t", index_col=None)
        print("    .sumstats.gz loaded, %d SNPs" % df_sumstats.shape[0])
        print("    " + gdreg.util.get_sys_info(sys_start_time))

    # Load --annot_file (lazy loading)
    if JOB in ["compute_score", "compute_avgr", "regress"]:
        print("# Loading --annot_file")
        dic_annot_path = {}
        dic_pannot_path = {}
        
        annot_file_list = []
        CHR0 = list(dic_data)[0]
        if ANNOT_FILE.endswith(".txt"):
            with open(ANNOT_FILE, "r") as f:
                for line in f:
                    line = line.strip()
                    if os.path.exists(line.replace("@", "%d" % CHR0)):
                        annot_file_list.append(line)
                    else:
                        print("    Skip: %s" % line)
        else:
            for line in ANNOT_FILE.split(','):
                line = line.strip()
                if os.path.exists(line.replace("@", "%d" % CHR0)):
                    annot_file_list.append(line)
                else:
                    print("    Skip: %s" % line)            

        for annot_file in annot_file_list:
            annot_file = annot_file.strip()
            if annot_file.endswith((".annot.gz", ".pannot_mat.npz")) is False:
                print("    Skip: %s" % annot_file)
                continue
            annot_name = gdreg.util.get_annot_name_from_file(annot_file)
            if annot_file.endswith(".annot.gz"):
                dic_annot_path[annot_name] = {}
                for CHR in dic_data:
                    fpath = annot_file.replace("@", "%d" % CHR)
                    if os.path.exists(fpath):
                        dic_annot_path[annot_name][CHR] = fpath
                CHR_set_annot = set(dic_annot_path[annot_name])
            if annot_file.endswith(".pannot_mat.npz"):
                dic_pannot_path[annot_name] = {}
                for CHR in dic_data:
                    fpath = annot_file.replace("@", "%d" % CHR)
                    if os.path.exists(fpath):
                        dic_pannot_path[annot_name][CHR] = fpath
                CHR_set_annot = set(dic_pannot_path[annot_name])

        # Check: all annots and pannots have the same set of CHRs
        for annot_name in dic_annot_path:
            err_msg = "Set of CHRs does not match for %s" % annot_name
            assert set(dic_annot_path[annot_name]) == CHR_set_annot, err_msg
        for annot_name in dic_pannot_path:
            err_msg = "Set of CHRs does not match for %s" % annot_name
            assert set(dic_pannot_path[annot_name]) == CHR_set_annot, err_msg
        print(
            "    Detected %d CHRs for all files: %s"
            % (len(CHR_set_annot), ",".join(["%d" % x for x in CHR_set_annot]))
        )

        # Check: annots have the same col_list across CHRs
        for annot_name in dic_annot_path:
            CHR0 = list(CHR_set_annot)[0]
            col_list = list(
                gdreg.util.read_annot(dic_annot_path[annot_name][CHR0], nrows=5)
            )
            for CHR in CHR_set_annot:
                temp_df = gdreg.util.read_annot(
                    dic_annot_path[annot_name][CHR], nrows=5
                )
                err_msg = "%s : columns mismatch between CHR%d and CHR%d" % (
                    annot_name,
                    CHR0,
                    CHR,
                )
                assert list(temp_df) == col_list, err_msg
            print("    %s: columms match for all CHRs. Containing" % annot_name)
            temp_str = ",".join([x for x in col_list if x.startswith("AN:")])
            print("        %s" % temp_str)

        # Check: pannots have the same shape as pvar file
        for annot_name in dic_pannot_path:
            CHR = np.random.choice(list(CHR_set_annot), size=1)[0]
            mat_G = gdreg.util.read_pannot_mat(dic_pannot_path[annot_name][CHR])
            err_msg = "(%s, CHR%d) : n_snp=%d, mismatch with --pgen_file" % (
                annot_name,
                CHR,
                mat_G.shape[0],
            )
            assert mat_G.shape[0] == dic_data[CHR]["pvar"].shape[0], err_msg
            print(
                "    %s (%d CHRs) : CHR%d dimension matches with .pvar"
                % (annot_name, len(dic_pannot_path[annot_name]), CHR)
            )

    # Load --avgr_file
    if JOB in ["regress"]:
        print("# Loading --avgr_file")
        if AVGR_FILE is not None:
            dic_avgr = {}
            for fpath in AVGR_FILE.split(","):
                fpath = fpath.strip()
                if os.path.exists(fpath):
                    temp_df = pd.read_csv(fpath, sep="\t", header=None)
                    dic_avgr.update(
                        {
                            x: y
                            for x, y in zip(temp_df[0], temp_df[1])
                            if x in dic_pannot_path
                        }
                    )
                else:
                    print("    Skip: %s" % fpath)
            print(
                "    Detected AVGR for %d/%d pannots"
                % (len(dic_avgr), len(dic_pannot_path))
            )
            for pAN in dic_pannot_path:
                if pAN not in dic_avgr:
                    print("    Missing AVGR for %s, pad 0" % pAN)
                    dic_avgr[pAN] = 0
        else:
            print("--avgr_file missing, add 0 to pannots")
            dic_avgr = {x: 0 for x in dic_pannot_path}
        for pAN in dic_avgr:
            print("    %s: %0.4f" % (pAN, dic_avgr[pAN]))

    # Load --score_file
    if JOB in ["regress"]:
        print("# Loading --score_file")
        # CHR_list_score containing all score files
        CHR_list_score = set(dic_data)
        for score_file in SCORE_FILE.split(","):
            score_file = score_file.strip()
            temp_list = []
            for CHR in CHR_list_score:
                if os.path.exists(score_file.replace("@", "%d" % CHR)):
                    temp_list.append(CHR)
            CHR_list_score = CHR_list_score & set(temp_list)
        CHR_list_score = sorted(CHR_list_score)
        print(
            "    Detected all score files for %d CHRs: %s"
            % (len(CHR_list_score), ",".join(["%d" % x for x in CHR_list_score]))
        )

        # Preload to get info
        df_score = None
        for score_file in SCORE_FILE.split(","):
            CHR0 = CHR_list_score[0]
            score_file = score_file.strip().replace("@", "%d" % CHR0)
            temp_df = pd.read_csv(score_file, sep="\t", index_col=None, nrows=100)
            temp_df.index = temp_df["SNP"]
            if df_score is None:
                df_score = temp_df.copy()
            else:
                col_list = [x for x in temp_df if x not in df_score]
                df_score = df_score.join(temp_df[col_list])

        AN_list_score = [x.replace("LD:", "") for x in df_score if x.startswith("LD:")]
        pAN_list_score = [
            x.replace("DLD:", "") for x in df_score if x.startswith("DLD:")
        ]

        print(
            "    Detected %d LD scores: %s"
            % (len(AN_list_score), ",".join(AN_list_score))
        )
        print(
            "    Detected %d DLD scores: %s"
            % (len(pAN_list_score), ",".join(pAN_list_score))
        )

        # Check: if all scores have annots/pannots
        AN_list, CHR0 = [], list(CHR_set_annot)[0]
        for annot_name in dic_annot_path:
            temp_df = gdreg.util.read_annot(dic_annot_path[annot_name][CHR0], nrows=5)
            AN_list.extend([x for x in temp_df if x.startswith("AN:")])
        pAN_list = list(dic_pannot_path)

        #         temp_list = list(set(AN_list_score) - set(AN_list))
        #         err_msg = "%d scores without annots: %s" % (
        #             len(temp_list),
        #             ",".join(temp_list),
        #         )
        #         assert len(temp_list) == 0, err_msg
        #         temp_list = list(set(pAN_list_score) - set(pAN_list))
        #         err_msg = "%d scores without pannots: %s" % (
        #             len(temp_list),
        #             ",".join(temp_list),
        #         )
        #         assert len(temp_list) == 0, err_msg

        drop_list_LD = ["LD:%s" % x for x in set(AN_list_score) - set(AN_list)]
        drop_list_DLD = ["DLD:%s" % x for x in set(pAN_list_score) - set(pAN_list)]
        if len(drop_list_LD) > 0:
            print(
                "    Drop %d LD scores without annots: %s"
                % (len(drop_list_LD), ",".join(drop_list_LD))
            )
        if len(drop_list_DLD) > 0:
            print(
                "    Drop %d DLD scores without pannots: %s"
                % (len(drop_list_DLD), ",".join(drop_list_DLD))
            )

        # Finally, load all score files
        df_score = None
        for score_file in SCORE_FILE.split(","):
            df_list = []
            for CHR in CHR_list_score:
                fpath = score_file.replace("@", "%d" % CHR)
                if os.path.exists(fpath):
                    temp_df = pd.read_csv(fpath, sep="\t", index_col=None)
                    col_list = [x for x in temp_df if x.startswith(("E", "LD", "DLD"))]
                    temp_df[col_list] = temp_df[col_list].astype(np.float32)
                    df_list.append(temp_df.copy())

            temp_df = pd.concat(df_list, axis=0)
            temp_df.index = temp_df["SNP"]
            if df_score is None:
                df_score = temp_df.copy()
            else:
                col_list = [x for x in temp_df if x not in df_score]
                df_score = df_score.join(temp_df[col_list])
            del temp_df

        df_score = df_score.loc[df_score.isna().sum(axis=1) == 0].copy()
        df_score.drop(columns=drop_list_LD + drop_list_DLD, inplace=True)
        df_score.sort_values(["CHR", "BP"], inplace=True)
        df_score.index = df_score["SNP"]
        n_snp = df_score.shape[0]
        n_LD = len([x for x in df_score if x.startswith("LD:")])
        n_DLD = len([x for x in df_score if x.startswith("DLD:")])
        print("    Loaded: %d SNPs, %d LD scores, %d DLD scores" % (n_snp, n_LD, n_DLD))
        print("    " + gdreg.util.get_sys_info(sys_start_time))

    ###########################################################################################
    ######                                  Computation                                  ######
    ###########################################################################################
    if JOB == "get_snp_block":
        print("# Running --job get_snp_block")
        block_size = 10000
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
                dic_data,
                pos_tar,
                pos_ref,
                verbose=True,
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
                ind_s_ref = np.searchsorted(v_bp, v_bp[ind_s] - 5.01e6, side="left")
                ind_s_ref = max(0, ind_s_ref - 1)
                ind_e_ref = np.searchsorted(
                    v_bp, v_bp[ind_e - 1] + 5.01e6, side="right"
                )
                ind_e_ref = min(n_snp, ind_e_ref + 1)

                pos_tar = [CHR, ind_s, ind_e]
                pos_ref = [CHR, ind_s_ref, ind_e_ref]
                mat_ld = gdreg.score.compute_ld(
                    dic_data,
                    pos_tar,
                    pos_ref,
                    verbose=True,
                )
                temp_mat = np.zeros([n_snp, ind_e - ind_s], dtype=np.float32)
                temp_mat[ind_s_ref:ind_e_ref, :] = mat_ld
                mat_ld_list.append(sp.sparse.csc_matrix(temp_mat))
                del temp_mat

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
            dic_annot_path=dic_annot_path,
            dic_pannot_path=dic_pannot_path,
            flag_cross_term=FLAG_CROSS_TERM,
            verbose=True,
            win_size=1e7,
            snp_range=snp_range,
        )

        col_list = [x for x in df_score if x.startswith(("E", "LD", "DLD"))]
        df_score[col_list] = df_score[col_list].astype(np.float32)

        df_score.to_csv(
            PREFIX_OUT + ".c%d_s%d_e%d_score.tsv.gz" % (CHR, START, END),
            sep="\t",
            index=False,
            compression="gzip",
        )
        print("    " + gdreg.util.get_sys_info(sys_start_time))

    if JOB == "compute_avgr":
        print("# Running --job compute_avgr")
        dic_avgr = gdreg.score.compute_avgr(dic_pannot_path, dic_ld_path, verbose=True)
        with open(PREFIX_OUT + ".avgr", "w") as f:
            for pAN in dic_avgr:
                f.write("%s\t%0.6f\n" % (pAN, dic_avgr[pAN]))
        print("    " + gdreg.util.get_sys_info(sys_start_time))

    if JOB == "regress":
        print("# Running --job regress")

        dic_res = gdreg.regress.estimate(
            dic_data,
            df_score,
            df_sumstats,
            dic_annot_path=dic_annot_path,
            dic_pannot_path=dic_pannot_path,
            dic_avgr=dic_avgr,
            flag_cross_term=FLAG_CROSS_TERM,
            flag_nofil_snp=FLAG_NOFIL_SNP,
            n_jn_block=100,
            verbose=True,
        )

        # Store the entire file and a summary df
        dbfile = open(PREFIX_OUT + ".pickle", "wb")
        pickle.dump(dic_res, dbfile)
        dbfile.close()
        dic_res["summary"]["tau"].to_csv(PREFIX_OUT + ".tau.tsv", sep="\t", index=False)
        dic_res["summary"]["rho"].to_csv(PREFIX_OUT + ".rho.tsv", sep="\t", index=False)

        print("    " + gdreg.util.get_sys_info(sys_start_time))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="gdreg")

    parser.add_argument("--job", type=str, required=True, help="One of [compute_ld]")
    parser.add_argument("--pgen_file", type=str, required=False, default=None)
    parser.add_argument(
        "--ld_file", type=str, required=False, default=None, help=".<range>_ld.npz"
    )
    parser.add_argument("--annot_file", type=str, required=False, default=None, help="Comma-separated file paths or .txt file with one line per file path.")
    parser.add_argument("--score_file", type=str, required=False, default=None)
    parser.add_argument("--sumstats_file", type=str, required=False, default=None)
    parser.add_argument("--avgr_file", type=str, required=False, default=None)
    parser.add_argument("--prefix_out", type=str, required=True)
    parser.add_argument(
        "--snp_range",
        type=str,
        default=None,
        help="c1_s20_e1701_r1, '_rall' for all ref CHRs",
    )
    parser.add_argument("--flag_full_ld", type=bool, default=False)
    parser.add_argument("--flag_cross_term", type=bool, default=False)
    parser.add_argument("--flag_nofil_snp", type=bool, default=False)

    args = parser.parse_args()
    main(args)
