import pandas as pd
import numpy as np
import time
import os
import argparse
import gdreg


"""
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
    EFF_FILE = args.eff_file
    PHEN_FILE = args.phen_file
    PREFIX_OUT = args.prefix_out
    MEMORY = args.memory
    RANDOM_SEED = args.random_seed

    # Parse and check arguments
    LEGAL_JOB_LIST = ["compute_phen", "compute_sumstats"]
    err_msg = "# run_gdreg: --job=%s not supported" % JOB
    assert JOB in LEGAL_JOB_LIST, err_msg

    if (EFF_FILE is None) & (JOB in ["compute_phen"]):
        raise ValueError("# run_simulation.py: --eff_file required for --job=%s" % JOB)
    if (PHEN_FILE is None) & (JOB in ["compute_sumstats"]):
        raise ValueError("# run_simulation.py: --phen_file required for --job=%s" % JOB)

    # Print input options
    header = gdreg.util.get_cli_head()
    header += "Call: run_simulation.py \\\n"
    header += "--job %s\\\n" % JOB
    header += "--pgen_file %s\\\n" % PGEN_FILE
    header += "--eff_file %s\\\n" % EFF_FILE
    header += "--phen_file %s\\\n" % PHEN_FILE
    header += "--prefix_out %s\\\n" % PREFIX_OUT
    header += "--memory %d\\\n" % MEMORY
    header += "--random_seed %d\n" % RANDOM_SEED
    print(header)

    ###########################################################################################
    ######                                   Data Loading                                ######
    ###########################################################################################
    # Load genotype data
    if JOB in ["compute_phen", "compute_sumstats"]:
        print("# Loading --pgen_file")
        dic_data = {}
        for CHR in range(1, 23):
            if os.path.exists(PGEN_FILE.replace("@", "%s" % CHR) + ".pgen"):
                dic_data[CHR] = gdreg.util.read_pgen(PGEN_FILE.replace("@", "%s" % CHR))

        print("    Genotype data for %d CHRs:" % len(dic_data))
        for CHR in dic_data:
            n_sample = dic_data[CHR]["psam"].shape[0]
            n_snp = dic_data[CHR]["pvar"].shape[0]
            print("        CHR%2d (%d samples %d SNPs)" % (CHR, n_sample, n_snp))
        print("    " + gdreg.util.get_sys_info(sys_start_time))

    # Load EFF_FILE
    if JOB in ["compute_phen"]:
        print("# Loading --eff_file")
        df_effect = pd.read_csv(EFF_FILE, sep="\t", index_col=None)

        n_snp = df_effect.shape[0]
        n_CHR = len(set(df_effect["CHR"]))
        h2g = (df_effect["EFF"] ** 2).sum()
        print("    %s SNPs from %d CHRs, h2g=%0.3f" % (n_snp, n_CHR, h2g))
        print("    " + gdreg.util.get_sys_info(sys_start_time))

    if JOB in ["compute_sumstats"]:
        print("# Loading --phen_file")
        df_phen = pd.read_csv(PHEN_FILE, sep="\t", index_col=None)
        phen_name = df_phen.columns[2]

        print("    %d samples, phen_name=%s" % (df_phen.shape[0], phen_name))
        print("    " + gdreg.util.get_sys_info(sys_start_time))

    ###########################################################################################
    ######                                  Computation                                  ######
    ###########################################################################################
    if JOB == "compute_phen":
        print("# Running --job compute_phen")
        df_phen = gdreg.simulate.simulate_phen(
            df_effect, dic_data, block_size=500, random_seed=RANDOM_SEED, verbose=True
        )
        df_phen.to_csv(PREFIX_OUT + ".phen", sep="\t", index=False)
        print("    " + gdreg.util.get_sys_info(sys_start_time))

    if JOB == "compute_sumstats":
        print("# Running --job compute_sumstats")
        df_sumstats = gdreg.simulate.compute_sumstats(
            df_phen, dic_data, block_size=500, verbose=True
        )
        df_sumstats.to_csv(PREFIX_OUT + ".sumstats.gz", sep="\t", index=False)
        print("    " + gdreg.util.get_sys_info(sys_start_time))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="gdreg")

    parser.add_argument(
        "--job",
        type=str,
        required=True,
        help="One of [compute_phen, compute_sumstats]",
    )
    parser.add_argument("--pgen_file", type=str, required=True)
    parser.add_argument("--eff_file", type=str, default=None)
    parser.add_argument("--phen_file", type=str, default=None)
    parser.add_argument("--prefix_out", type=str, required=True)
    parser.add_argument("--memory", type=str, default=128)
    parser.add_argument("--random_seed", type=int, default=0)

    args = parser.parse_args()
    main(args)
