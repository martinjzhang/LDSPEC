import gzip
import pandas as pd
import numpy as np
import time
import argparse
import gdreg

def main(args):
    
    sys_start_time = time.time()
    
    SNP_RANGE_FILE = "/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/ukb_imp_v3.snp_range.txt"

    PGEN_FILE = args.pgen_file
    LB = args.lb
    UB = args.ub
    LD_FILE = args.ld_file
    OUT_PATH = args.out_path

    print("ldp5_proxy_10000 .pannot_hr.gz")
    print("--pgen_file %s" % PGEN_FILE)
    print("--lb %s" % LB)
    print("--ub %s" % UB)
    print("--ld_file %s" % LD_FILE)
    print("--out_path %s" % OUT_PATH)

    df_snp_chr = gdreg.util.read_pgen(PGEN_FILE)["pvar"]
    df_snp_chr = df_snp_chr[["CHR", "SNP", "BP"]].copy()
    CHR = df_snp_chr["CHR"][0]
    df_snp_chr['MAF'] = gdreg.util.read_pgen(PGEN_FILE)["afreq"]["MAF"].astype(np.float32)
    dic_bp = {x:y for x,y in zip(df_snp_chr["SNP"], df_snp_chr["BP"])}

    # SNP pairs with LD>0.5 and dist<1e6
    snp_range_list = list(pd.read_csv(SNP_RANGE_FILE, header=None)[0])
    snp_range_list_chr = [x for x in snp_range_list if x.startswith("c%d_" % CHR)]

    snp_list = list(df_snp_chr["SNP"])
    snp_list1 = []
    snp_list2 = []

    for snp_range in snp_range_list_chr:
        print("    snp_range=%s" % snp_range)
        mat_ld, dic_range = gdreg.util.read_ld(LD_FILE.replace("@", snp_range))

        for i in range(dic_range["start"], dic_range["end"]):
            v_ld = mat_ld[:, i-dic_range["start"]].toarray()
            for j in range(i+1, len(snp_list)):
                ld = v_ld[j]
                dist = dic_bp[snp_list[j]] - dic_bp[snp_list[i]] 
                if (ld>=0.5) & (dist>=LB) & (dist<UB):
                    snp_list1.append(snp_list[i])
                    snp_list2.append(snp_list[j])
                elif dist>=UB:
                    break
                    
            if i%1000 == 0:
                print("        %d/%d SNPs" % (i-dic_range["start"], mat_ld.shape[1]))

    dic_maf_snp = {
        "common" : set(df_snp_chr["SNP"][df_snp_chr["MAF"].values >= 0.05]),
        "lf" : set(df_snp_chr["SNP"][(df_snp_chr["MAF"].values >= 0.005) & (df_snp_chr["MAF"].values < 0.05)]),
        "rare" : set(df_snp_chr["SNP"][df_snp_chr["MAF"].values <= 0.005]),
    }

    for term1,term2 in [
        ['common', 'common'], ['common', 'lf'], ['lf', 'lf'],
    ]:
        temp_snp_list1 = []
        temp_snp_list2 = []
        for snp1,snp2 in zip(snp_list1, snp_list2):
            flag1 = (snp1 in dic_maf_snp[term1]) & (snp2 in dic_maf_snp[term2])
            flag2 = (snp1 in dic_maf_snp[term2]) & (snp2 in dic_maf_snp[term1])
            if flag1 | flag2:
                temp_snp_list1.append(snp1)
                temp_snp_list2.append(snp2)   

        snp_pair_list = [(x,y) for x,y in zip(temp_snp_list1, temp_snp_list2)]
        if len(snp_pair_list) > 10:
            gdreg.util.write_pannot_mat(
                snp_pair_list, list(df_snp_chr["SNP"]), OUT_PATH + "/ldp5_proxy_%d_%d_%s_%s.chr%d" % (LB, UB, term1, term2, CHR)
            )
        print('pAN:ldp5_proxy_%d_%d_%s_%s' % (LB, UB, term1, term2), 'size=%d'% len(temp_snp_list1)) 

    print('# Finished, time=%0.1fs'%(time.time() - sys_start_time))
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='gdreg')
    parser.add_argument('--pgen_file', type=str, default=None)
    parser.add_argument('--lb', type=int, default=None)
    parser.add_argument('--ub', type=int, default=None)
    parser.add_argument('--ld_file', type=str, default=None)
    parser.add_argument('--out_path', type=str, default=None)
    
    args = parser.parse_args()
    main(args)