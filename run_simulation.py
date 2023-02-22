import pandas as pd
import numpy as np
import time
import os
import argparse
import gdreg


"""
Job description
----------------

compute_eff : simulate and summerize SNP effects 
    - Input : --pgen_file | --config_file | --annot_file | --prefix_out
    - Output : list of snp ranges (in the format of "snp_range")

compute_phen : compute .phen file
    - Input : --job | --pgen_file | --config_file | --annot_file |  [--random_seed] | [--flag_full_ld]
    - Output : LD matrix between a set of SNPs and all other SNPs on the same chromosome.
    
compute_sumstats : compute .sumstats file
    - Input : --job | --pgen_file | --ld_file | --annot_file | --prefix_out | [--random_seed] 
    | [--flag_cross_term]
    - Output : LD and DLD scores.
    
TODO
----
- Add truth of "cov(ac;proxy)" and "cor(ac;proxy)"
- 
"""


def main(args):
    sys_start_time = time.time()

    ###########################################################################################
    ######                                    Parse Options                              ######
    ###########################################################################################

    JOB = args.job
    PGEN_FILE = args.pgen_file
    CONFIG_FILE = args.config_file
    ANNOT_FILE = args.annot_file
    EFF_FILE = args.eff_file
    PHEN_FILE = args.phen_file
    PREFIX_OUT = args.prefix_out
    RANDOM_SEED = int(args.random_seed)
    FLAG_BW_SPARSE = bool(args.flag_bw_sparse)

    # Parse and check arguments
    LEGAL_JOB_LIST = [
        "simulate",
        "compute_phen",
        "compute_sumstats",
    ]
    err_msg = "# run_gdreg: --job=%s not supported" % JOB
    assert JOB in LEGAL_JOB_LIST, err_msg

    if (PHEN_FILE is None) & (JOB in ["compute_sumstats"]):
        raise ValueError("# run_simulation.py: --phen_file required for --job=%s" % JOB)

    # Print input options
    header = gdreg.util.get_cli_head()
    header += "Call: run_simulation.py \\\n"
    header += "--job %s\\\n" % JOB
    header += "--pgen_file %s\\\n" % PGEN_FILE
    header += "--config_file %s\\\n" % CONFIG_FILE
    header += "--annot_file %s\\\n" % ANNOT_FILE
    header += "--eff_file %s\\\n" % EFF_FILE
    header += "--phen_file %s\\\n" % PHEN_FILE
    header += "--prefix_out %s\\\n" % PREFIX_OUT
    header += "--random_seed %d\\\n" % RANDOM_SEED
    header += "--flag_bw_sparse %s\n" % FLAG_BW_SPARSE
    print(header)

    ###########################################################################################
    ######                                   Data Loading                                ######
    ###########################################################################################
    # Load genotype data
    if JOB in ["simulate", "compute_phen", "compute_sumstats"]:
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

    # Load --annot_file (lazy loading)
    if JOB in ["simulate", "compute_phen"]:
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
            for line in ANNOT_FILE.split(","):
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

    # Load config
    if JOB in ["simulate", "compute_phen"]:
        print("# Loading --config_file")
        temp_df = pd.read_csv(CONFIG_FILE, sep="\t", header=None)
        dic_config = {x: y for x, y in zip(temp_df[0], temp_df[1])}
        for col in ["h2g", "p_causal", "alpha"]:
            assert col in dic_config, "%s not in --config_file"
            print("    %s=%0.3f" % (col, dic_config[col]))

        AN_list, CHR0 = [], list(CHR_set_annot)[0]
        for annot_name in dic_annot_path:
            temp_df = gdreg.util.read_annot(dic_annot_path[annot_name][CHR0], nrows=5)
            AN_list.extend([x for x in temp_df if x.startswith("AN:")])
        pAN_list = list(dic_pannot_path)

        dic_coef = {
            x: dic_config[x]
            for x in dic_config
            if x not in ["h2g", "p_causal", "alpha"]
        }
        for annot in dic_coef:
            if annot not in ["h2g", "p_causal", "alpha"] + AN_list + pAN_list:
                print("    %s not in --annot_file" % annot)

        print("    %s" % ", ".join(["%s (%0.2f)" % (x, dic_coef[x]) for x in dic_coef]))

    # Load EFF_FILE
    if JOB in ["compute_phen"]:
        print("# Loading --eff_file")
        df_effect = pd.read_csv(EFF_FILE, sep="\t", index_col=None)
        print(
            "    %d SNPs, h2=%0.3f"
            % (df_effect.shape[0], (df_effect["EFF"] ** 2).sum())
        )
        print("    " + gdreg.util.get_sys_info(sys_start_time))

    # Load PHEN_FILE
    if JOB in ["compute_sumstats"]:
        print("# Loading --phen_file")
        df_phen = pd.read_csv(PHEN_FILE, sep="\t", index_col=None)
        phen_name = df_phen.columns[2]
        print("    %d samples, phen_name=%s" % (df_phen.shape[0], phen_name))
        print("    " + gdreg.util.get_sys_info(sys_start_time))

    ###########################################################################################
    ######                                  Computation                                  ######
    ###########################################################################################
    if JOB == "simulate":
        # Simulate SNP effects
        df_effect = gdreg.simulate.simulate_snp_effect(
            dic_data,
            dic_coef,
            dic_annot_path=dic_annot_path,
            dic_pannot_path=dic_pannot_path,
            h2g=dic_config["h2g"],
            alpha=dic_config["alpha"],
            p_causal=dic_config["p_causal"],
            block_size=100,
            flag_bw_sparse=FLAG_BW_SPARSE,
            random_seed=RANDOM_SEED,
            verbose=True,
        )
        df_effect.to_csv(PREFIX_OUT + ".eff.gz", sep="\t", index=False)
        print("    " + gdreg.util.get_sys_info(sys_start_time))

        # Compute .phen
        df_effect_ = df_effect.copy()
        df_phen = gdreg.simulate.simulate_phen(
            dic_data,
            dic_coef,
            df_effect_,
            dic_annot_path=dic_annot_path,
            block_size=500,
            random_seed=RANDOM_SEED + 42,
            verbose=True,
        )
        df_phen.to_csv(PREFIX_OUT + ".phen", sep="\t", index=False)
        print("    " + gdreg.util.get_sys_info(sys_start_time))

        # Summarize SNP effects
        df_effect_ = df_effect.copy()
        df_phen_ = df_phen.copy()
        df_sum_tau, df_sum_rho = gdreg.simulate.summarize_snp_effect(
            dic_data,
            dic_coef,
            df_effect_,
            df_phen_,
            dic_annot_path=dic_annot_path,
            dic_pannot_path=dic_pannot_path,
            block_size=1000,
            verbose=True,
        )
        df_sum_tau.to_csv(PREFIX_OUT + ".eff_tau.tsv", sep="\t", index=False)
        df_sum_rho.to_csv(PREFIX_OUT + ".eff_rho.tsv", sep="\t", index=False)
        print("    " + gdreg.util.get_sys_info(sys_start_time))

        # Compute .sumstats
        df_phen_ = df_phen.copy()
        df_sumstats = gdreg.simulate.compute_sumstats(
            df_phen_, dic_data, block_size=500, verbose=True
        )
        df_sumstats.to_csv(PREFIX_OUT + ".sumstats.gz", sep="\t", index=False)
        print("    " + gdreg.util.get_sys_info(sys_start_time))

    if JOB == "compute_phen":
        df_phen = gdreg.simulate.simulate_phen(
            dic_data,
            dic_coef,
            df_effect,
            dic_annot_path=dic_annot_path,
            block_size=500,
            random_seed=RANDOM_SEED + 42,
            verbose=True,
        )
        df_phen.to_csv(PREFIX_OUT + ".phen", sep="\t", index=False)
        print("    " + gdreg.util.get_sys_info(sys_start_time))

    if JOB == "compute_sumstats":
        df_sumstats = gdreg.simulate.compute_sumstats(
            df_phen, dic_data, block_size=500, verbose=True
        )
        df_sumstats.to_csv(PREFIX_OUT + ".sumstats.gz", sep="\t", index=False)
        print("    " + gdreg.util.get_sys_info(sys_start_time))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="gdreg")

    parser.add_argument("--job", type=str, default="simulate")
    parser.add_argument("--pgen_file", type=str, required=True)
    parser.add_argument("--config_file", type=str, default=None)
    parser.add_argument("--annot_file", type=str, default=None)
    parser.add_argument("--eff_file", type=str, default=None)
    parser.add_argument("--phen_file", type=str, default=None)
    parser.add_argument("--prefix_out", type=str, required=True)
    parser.add_argument("--flag_bw_sparse", type=bool, default=False)
    parser.add_argument("--random_seed", type=int, default=0)

    args = parser.parse_args()
    main(args)
