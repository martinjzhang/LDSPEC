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

compute_ld : compute LD matrix.
    - Input : --job | --pgen_file | --prefix_out | --snp_range | [--memory] ｜ [--random_seed]
        | [--flag_full_ld]
    - Output : LD matrix between a set of SNPs and all other SNPs on the same chromosome.
    
regress : infer parameters \tau and \rho.
    - Input : --job | --pgen_file | --ld_file | --sumstats_file | --annot_file | --prefix_out 
        | [--memory]
    - Output : GDREG result.
    
TODO
----
- compute_ld : provide option of computing LD for only a range of genome (used when data is large).
- load dic_ld only when using it

"""


def main(args):
    sys_start_time = time.time()

    ###########################################################################################
    ######                                    Parse Options                              ######
    ###########################################################################################

    JOB = args.job
    PGEN_FILE = args.pgen_file
    LD_FILE = args.ld_file
    SUMSTATS_FILE = args.sumstats_file
    ANNOT_FILE = args.annot_file
    PREFIX_OUT = args.prefix_out
    MEMORY = args.memory
    RANDOM_SEED = args.random_seed
    SNP_RANGE = args.snp_range
    FLAG_FULL_LD = args.flag_full_ld

    # Parse and check arguments
    LEGAL_JOB_LIST = ["compute_ld", "compute_score", "regress"]
    err_msg = "# run_gdreg: --job=%s not supported" % JOB
    assert JOB in LEGAL_JOB_LIST, err_msg

    if JOB in ["compute_score", "regress"]:
        assert LD_FILE is not None, "--ld_file required for --job=%s" % JOB
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
    if JOB in ["compute_ld", "compute_score", "regress"]:
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
    if JOB in ["compute_score", "regress"]:
        print("# Loading --ld_file")
        dic_ld = {}
        for CHR in dic_data:
            err_msg = "--ld_file missing for CHR%d" % CHR
            assert os.path.exists(LD_FILE.replace("@", "%s" % CHR)), err_msg
            if LD_FILE.endswith(".full_ld.npy"):
                dic_ld[CHR] = np.load(LD_FILE.replace("@", "%s" % CHR))
            elif LD_FILE.endswith(".ld.npz"):
                dic_ld[CHR] = sp.sparse.load_npz(LD_FILE.replace("@", "%s" % CHR))
            err_msg = "CHR%2d n_snp=%d, mismatch with --pgen_file" % (
                CHR,
                dic_ld[CHR].shape[0],
            )
            assert dic_ld[CHR].shape[0] == dic_data[CHR]["pvar"].shape[0], err_msg
        print("    LD info loaded, matching --pgen_file")
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
    if JOB == "compute_ld":
        print("# Running --job compute_ld")
        if FLAG_FULL_LD:
            CHR, CHR_REF = DIC_RANGE["chr"], DIC_RANGE["chr_ref"][0]
            pos_tar = [CHR, 0, dic_data[CHR]["pvar"].shape[0]]
            pos_ref = [CHR_REF, 0, dic_data[CHR_REF]["pvar"].shape[0]]
            mat_ld = gdreg.score.compute_ld(
                dic_data, pos_tar, pos_ref, verbose=True, memory=MEMORY
            )
            np.save(PREFIX_OUT + ".full_ld", mat_ld)
        else:
            CHR = DIC_RANGE["chr"]
            n_snp = dic_data[CHR]["pvar"].shape[0]
            v_bp = dic_data[CHR]["pvar"]["BP"].values

            block_size = 1000
            n_block = np.ceil(n_snp / block_size).astype(int)
            mat_ld_list = []
            for i_block in range(n_block):
                print("block %d/%d %s"%(i_block, n_block, gdreg.util.get_sys_info(sys_start_time)))
                ind_s = i_block * block_size 
                ind_e = min((i_block + 1) * block_size, n_snp)        
                ind_s_ref = np.searchsorted(v_bp, v_bp[ind_s] - 5.1e6, side="left")
                ind_s_ref = max(0, ind_s_ref-5)
                ind_e_ref = np.searchsorted(v_bp, v_bp[ind_e - 1] + 5.1e6, side="right")
                ind_e_ref = min(n_snp, ind_e_ref+5)

                pos_tar = [CHR, ind_s, ind_e]
                pos_ref = [CHR, ind_s_ref, ind_e_ref]
                mat_ld = gdreg.score.compute_ld(
                    dic_data, pos_tar, pos_ref, verbose=True, memory=MEMORY, 
                )
                temp_mat = np.zeros([n_snp, ind_e-ind_s], dtype=np.float32)
                temp_mat[ind_s_ref:ind_e_ref, :] = mat_ld
                mat_ld_list.append(sp.sparse.csc_matrix(temp_mat))

            mat_ld = sp.sparse.hstack(mat_ld_list, format='csc')
            sp.sparse.save_npz(PREFIX_OUT + ".ld", mat_ld)
        print("    " + gdreg.util.get_sys_info(sys_start_time))
            

#     if JOB == "compute_score":
#         pass

    if JOB == "regress":
        print("# Running --job regress")

        dic_res = gdreg.regress.estimate(
            dic_data,
            df_sumstats,
            dic_ld,
            df_annot,
            pannot_list=pannot_list,
            pannot_hr_list=pannot_hr_list,
            n_jn_block=100,
            sym_non_pAN="non-pAN",
            win_size=int(1e7),
            memory=MEMORY,
            verbose=True,
            n_iter=10,
        )

        # Store the entire file and a summary df
        dbfile = open(PREFIX_OUT+'.pickle', 'wb')      
        pickle.dump(dic_res, dbfile)                     
        dbfile.close()
        for res in dic_res:
            dic_res[res]['summary'].to_csv(PREFIX_OUT+'_res%s.tsv' % res, sep='\t', index=False)
        dic_res[res]['summary'].to_csv(PREFIX_OUT+'.tsv', sep='\t', index=False)
        print("    " + gdreg.util.get_sys_info(sys_start_time))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="gdreg")

    parser.add_argument("--job", type=str, required=True, help="One of [compute_ld]")
    parser.add_argument("--pgen_file", type=str, required=True)
    parser.add_argument("--ld_file", type=str, required=False, default=None)
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
