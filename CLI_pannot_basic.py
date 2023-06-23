import gzip
import pandas as pd
import numpy as np
import time
import argparse
import ldspec

"""
Create proximity-based pannot annotations.

Basic proximity-based annotation: 
    --pgen_file 
    --dist_lb
    --dist_ub
    --out_path
    --maf_ratio_thres

Stratifying by LD
    --ld_file : LDSPEC-computed LD matrices
    --snp_range_file : corresponding LDSEPC SNP-range file
    --ld_lb
    --ld_ub

Stratifying by MAF:
    --maf_bin_file : 
        Example: 
        common 0.05 0.5 
        lf 0.005 0.05
"""


def main(args):

    sys_start_time = time.time()

    PGEN_FILE = args.pgen_file
    LD_FILE = args.ld_file
    SNP_RANGE_FILE = args.snp_range_file
    DIST_LB = int(args.dist_lb)
    DIST_UB = int(args.dist_ub)
    LD_LB = float(args.ld_lb)
    LD_UB = float(args.ld_ub)
    MAF_BIN_FILE = args.maf_bin_file
    MAF_RATIO_THRES = int(args.maf_ratio_thres)
    FLAG_MAF_BLOCK = args.flag_maf_block
    OUT_PATH = args.out_path

    print("--pgen_file %s" % PGEN_FILE)
    print("--ld_file %s" % LD_FILE)
    print("--snp_range_file %s" % SNP_RANGE_FILE)
    print("--dist_lb %s" % DIST_LB)
    print("--dist_ub %s" % DIST_UB)
    print("--ld_lb %s" % LD_LB)
    print("--ld_ub %s" % LD_UB)
    print("--maf_bin_file %s" % MAF_BIN_FILE)
    print("--maf_ratio_thres %s" % MAF_RATIO_THRES)
    print("--flag_maf_block %s" % FLAG_MAF_BLOCK)
    print("--out_path %s" % OUT_PATH)

    # Data loading
    df_snp_chr = ldspec.util.read_pgen(PGEN_FILE)["pvar"]
    df_snp_chr = df_snp_chr[["CHR", "SNP", "BP"]].copy()
    CHR = df_snp_chr["CHR"][0]
    df_snp_chr["AF"] = ldspec.util.read_pgen(PGEN_FILE)["afreq"]["MAF"].astype(
        np.float32
    )
    df_snp_chr["MAF"] = [min(x, 1 - x) for x in df_snp_chr["AF"]]
    dic_bp = {x: y for x, y in zip(df_snp_chr["SNP"], df_snp_chr["BP"])}
    dic_maf = {x: y for x, y in zip(df_snp_chr["SNP"], df_snp_chr["MAF"])}
    print("n_snp", df_snp_chr.shape[0])

    # Data loading: stratifying by LD
    if LD_FILE is not None:
        snp_range_list = list(pd.read_csv(SNP_RANGE_FILE, header=None)[0])
        snp_range_list_chr = [x for x in snp_range_list if x.startswith("c%d_" % CHR)]

    # Data loading: stratifying by MAF
    if MAF_BIN_FILE is not None:
        df_mbin = pd.read_csv(MAF_BIN_FILE, delim_whitespace=True, header=None)
    else:
        df_mbin = pd.DataFrame(data={0: ["all"], 1: [0], 2: [1]})
    if FLAG_MAF_BLOCK:  # Disable MAF_RATIO_THRES
        MAF_RATIO_THRES = 1e6
    print(df_mbin)

    # Computation
    snp_list = list(df_snp_chr["SNP"])
    snp_list1 = []
    snp_list2 = []

    if LD_FILE is None:
        # No LD stratification
        for i in range(len(snp_list)):
            for j in range(i + 1, len(snp_list)):
                dist = dic_bp[snp_list[j]] - dic_bp[snp_list[i]]
                if dist > DIST_UB:
                    break
                maf_ratio = dic_maf[snp_list[i]] / dic_maf[snp_list[j]]
                maf_ratio = max(maf_ratio, 1 / maf_ratio)
                if (dist >= DIST_LB) & (dist < DIST_UB) & (maf_ratio < MAF_RATIO_THRES):
                    snp_list1.append(snp_list[i])
                    snp_list2.append(snp_list[j])
    else:
        # LD stratification
        for snp_range in snp_range_list_chr:
            print("    snp_range=%s" % snp_range)
            mat_ld, dic_range = ldspec.util.read_ld(LD_FILE.replace("@", snp_range))
            for i in range(dic_range["start"], dic_range["end"]):
                v_ld = mat_ld[:, i - dic_range["start"]].toarray().flatten()
                for j in range(i + 1, len(snp_list)):
                    dist = dic_bp[snp_list[j]] - dic_bp[snp_list[i]]
                    if dist > DIST_UB:
                        break
                    ld = v_ld[j]
                    maf_ratio = dic_maf[snp_list[i]] / dic_maf[snp_list[j]]
                    maf_ratio = max(maf_ratio, 1 / maf_ratio)
                    if (
                        (ld >= LD_LB)
                        & (ld <= LD_UB)
                        & (dist >= DIST_LB)
                        & (dist < DIST_UB)
                        & (maf_ratio < MAF_RATIO_THRES)
                    ):
                        snp_list1.append(snp_list[i])
                        snp_list2.append(snp_list[j])

    print("n_pair", len(snp_list1))

    # MAF-bin & write files
    for mbin, maf_lb, maf_ub in zip(df_mbin[0], df_mbin[1], df_mbin[2]):
        temp_snp_list1 = []
        temp_snp_list2 = []
        if FLAG_MAF_BLOCK:  # snps in the same MAF bin
            for snp1, snp2 in zip(snp_list1, snp_list2):
                if (
                    (dic_maf[snp1] >= maf_lb)
                    & (dic_maf[snp1] < maf_ub)
                    & (dic_maf[snp2] >= maf_lb)
                    & (dic_maf[snp2] < maf_ub)
                ):
                    temp_snp_list1.append(snp1)
                    temp_snp_list2.append(snp2)
        else:  # snps with geometric mean in the MAF bin
            for snp1, snp2 in zip(snp_list1, snp_list2):
                mean_maf = np.sqrt(dic_maf[snp1] * dic_maf[snp2])
                if (mean_maf >= maf_lb) & (mean_maf < maf_ub):
                    temp_snp_list1.append(snp1)
                    temp_snp_list2.append(snp2)

        snp_pair_list = [(x, y) for x, y in zip(temp_snp_list1, temp_snp_list2)]
        if len(snp_pair_list) > 10:
            str_prox = "proxy_%d_%d" % (DIST_LB, DIST_UB)            
            LDLB_str = (
                "n%d" % (int(-LD_LB * 100)) if LD_LB < 0 else "p%d" % (int(LD_LB * 100))
            )
            LDUB_str = (
                "n%d" % (int(-LD_UB * 100)) if LD_UB < 0 else "p%d" % (int(LD_UB * 100))
            )
            str_ld = (
                "ld_full" if LD_FILE is None else "ld_%s_%s" % (LDLB_str, LDUB_str)
            )
            str_maf = (
                "maf_%s_block" % mbin if FLAG_MAF_BLOCK else "maf_%s_geomean" % mbin
            )
            file_name = "%s.%s.%s.c%d" % (str_prox, str_ld, str_maf, CHR)
            print("%-50s" % file_name, "n_pair=%d" % len(snp_pair_list))
            ldspec.util.write_pannot_mat(
                snp_pair_list, list(df_snp_chr["SNP"]), OUT_PATH + "/" + file_name
            )

    print("# Finished, time=%0.1fs" % (time.time() - sys_start_time))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ldspec")
    parser.add_argument("--pgen_file", type=str, default=None, required=True)
    parser.add_argument("--ld_file", type=str, default=None, required=False)
    parser.add_argument("--snp_range_file", type=str, default=None, required=False)
    parser.add_argument("--dist_lb", type=str, default="0", required=False)
    parser.add_argument("--dist_ub", type=str, default="1000", required=False)
    parser.add_argument("--ld_lb", type=str, default="-1", required=False)
    parser.add_argument("--ld_ub", type=str, default="1", required=False)
    parser.add_argument("--maf_bin_file", type=str, default=None, required=False)
    parser.add_argument("--maf_ratio_thres", type=str, default="5", required=False)
    parser.add_argument("--flag_maf_block", type=bool, default=True, required=False)
    parser.add_argument("--out_path", type=str, default=None, required=False)

    args = parser.parse_args()
    main(args)
