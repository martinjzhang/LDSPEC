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
- Double check file consistency 
- 
"""


def main(args):
    sys_start_time = time.time()

    ###########################################################################################
    ######                                    Parse Options                              ######
    ###########################################################################################

    PGEN_FILE = args.pgen_file
    CONFIG_FILE = args.config_file
    ANNOT_FILE = args.annot_file
    PREFIX_OUT = args.prefix_out
    RANDOM_SEED = int(args.random_seed)

    # Print input options
    header = gdreg.util.get_cli_head()
    header += "Call: run_simulation.py \\\n"
    header += "--pgen_file %s\\\n" % PGEN_FILE
    header += "--config_file %s\\\n" % CONFIG_FILE
    header += "--annot_file %s\\\n" % ANNOT_FILE
    header += "--prefix_out %s\\\n" % PREFIX_OUT
    header += "--random_seed %d\n" % RANDOM_SEED
    print(header)

    ###########################################################################################
    ######                                   Data Loading                                ######
    ###########################################################################################
    # Load genotype data
    print("# Loading --pgen_file")
    dic_data = {}
    if "@" not in PGEN_FILE:  # Load single CHR
        temp_dic = gdreg.util.read_pgen(PGEN_FILE)
        dic_data[temp_dic["pvar"]["CHR"][0]] = temp_dic.copy()
    else:
        for CHR in range(1, 23):  # Check all 23 CHRs
            if os.path.exists(PGEN_FILE.replace("@", "%s" % CHR) + ".pgen"):
                dic_data[CHR] = gdreg.util.read_pgen(PGEN_FILE.replace("@", "%s" % CHR))

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
    print("# Loading --annot_file")
    dic_annot_path = {}
    dic_pannot_path = {}

    for annot_file in ANNOT_FILE.split(","):
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
            temp_df = gdreg.util.read_annot(dic_annot_path[annot_name][CHR], nrows=5)
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
        x: dic_config[x] for x in dic_config if x not in ["h2g", "p_causal", "alpha"]
    }
    for annot in dic_coef:
        if annot not in ["h2g", "p_causal", "alpha"] + AN_list + pAN_list:
            print("    %s not in --annot_file" % annot)

    print("    %s" % ", ".join(["%s (%0.2f)" % (x, dic_coef[x]) for x in dic_coef]))

    ###########################################################################################
    ######                                  Computation                                  ######
    ###########################################################################################
    # Simulate SNP effects
    df_effect = gdreg.simulate.simulate_snp_effect(
        dic_data,
        dic_coef,
        dic_annot_path=dic_annot_path,
        dic_pannot_path=dic_pannot_path,
        h2g=dic_config["h2g"],
        alpha=dic_config["alpha"],
        p_causal=dic_config["p_causal"],
        block_size=1000,
        random_seed=RANDOM_SEED,
        verbose=True,
    )
    df_effect.to_csv(PREFIX_OUT + ".eff.gz", sep="\t", index=False)
    print("    " + gdreg.util.get_sys_info(sys_start_time))

    # Compute .phen
    df_phen = gdreg.simulate.simulate_phen(
        dic_data,
        dic_coef,
        df_effect,
        dic_annot_path=dic_annot_path,
        block_size=500,
        random_seed=RANDOM_SEED,
        verbose=True,
    )
    df_phen.to_csv(PREFIX_OUT + ".phen", sep="\t", index=False)
    print("    " + gdreg.util.get_sys_info(sys_start_time))

    # Summarize SNP effects
    df_sum_tau, df_sum_rho = gdreg.simulate.summarize_snp_effect(
        dic_data,
        dic_coef,
        df_effect,
        df_phen,
        dic_annot_path=dic_annot_path,
        dic_pannot_path=dic_pannot_path,
        block_size=1000,
        verbose=True,
    )
    df_sum_tau.to_csv(PREFIX_OUT + ".eff_tau.tsv", sep="\t", index=False)
    df_sum_rho.to_csv(PREFIX_OUT + ".eff_rho.tsv", sep="\t", index=False)
    print("    " + gdreg.util.get_sys_info(sys_start_time))

    # Compute .sumstats
    df_sumstats = gdreg.simulate.compute_sumstats(
        df_phen, dic_data, block_size=500, verbose=True
    )
    df_sumstats.to_csv(PREFIX_OUT + ".sumstats.gz", sep="\t", index=False)
    print("    " + gdreg.util.get_sys_info(sys_start_time))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="gdreg")

    parser.add_argument("--pgen_file", type=str, required=True)
    parser.add_argument("--config_file", type=str, default=None)
    parser.add_argument("--annot_file", type=str, default=None)
    parser.add_argument("--prefix_out", type=str, required=True)
    parser.add_argument("--random_seed", type=int, default=0)

    args = parser.parse_args()
    main(args)
