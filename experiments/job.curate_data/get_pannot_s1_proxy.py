import gzip
import pandas as pd
import numpy as np
import time
import argparse
import gdreg

def main(args):
    
    sys_start_time = time.time()
    
    PGEN_FILE = args.pgen_file
    LB = args.lb
    UB = args.ub
    OUT_PATH = args.out_path
    
    print("proxy .pannot_hr.gz")
    print("--pgen_file %s" % PGEN_FILE)
    print("--lb %s" % LB)
    print("--ub %s" % UB)
    print("--out_path %s" % OUT_PATH)
    
    df_snp_chr = gdreg.util.read_pgen(PGEN_FILE)["pvar"]
    df_snp_chr = df_snp_chr[["CHR", "SNP", "BP"]].copy()
    CHR = df_snp_chr["CHR"][0]
    df_snp_chr['MAF'] = gdreg.util.read_pgen(PGEN_FILE)["afreq"]["MAF"].astype(np.float32)
    dic_bp = {x:y for x,y in zip(df_snp_chr["SNP"], df_snp_chr["BP"])}

    # Nearby SNP pairs
    snp_list = list(df_snp_chr["SNP"])
    snp_list1 = []
    snp_list2 = []
    for i in range(len(snp_list)):
        for j in range(i+1, len(snp_list)):
            dist = dic_bp[snp_list[j]] - dic_bp[snp_list[i]] 
            if (dist>=LB) & (dist<UB):
                snp_list1.append(snp_list[i])
                snp_list2.append(snp_list[j])
            elif dist>=UB:
                break

    dic_maf_snp = {
        "common" : set(df_snp_chr["SNP"][df_snp_chr["MAF"].values >= 0.05]),
        "lf" : set(df_snp_chr["SNP"][(df_snp_chr["MAF"].values >= 0.005) & (df_snp_chr["MAF"].values < 0.05)]),
        "rare" : set(df_snp_chr["SNP"][df_snp_chr["MAF"].values <= 0.005]),
    }

    for term1,term2 in [['common', 'common'], ['common', 'lf'], ['common', 'rare'], 
                 ['lf', 'lf'], ['lf', 'rare'], ['rare', 'rare']]:
        temp_snp_list1 = []
        temp_snp_list2 = []
        for snp1,snp2 in zip(snp_list1, snp_list2):
            flag1 = (snp1 in dic_maf_snp[term1]) & (snp2 in dic_maf_snp[term2])
            flag2 = (snp1 in dic_maf_snp[term2]) & (snp2 in dic_maf_snp[term1])
            if flag1 | flag2:
                temp_snp_list1.append(snp1)
                temp_snp_list2.append(snp2)   
                
        snp_pair_list = [(x,y) for x,y in zip(temp_snp_list1, temp_snp_list2)]
        gdreg.util.write_pannot_mat(
            snp_pair_list, list(df_snp_chr["SNP"]), OUT_PATH + "/proxy_%d_%d_%s_%s.chr%d" % (LB, UB, term1, term2, CHR)
        )
        print('pAN:gene_%s_%s' % (term1, term2), 'size=%d'% len(temp_snp_list1)) 

#         df_pannot = pd.DataFrame(data={
#             'CHR' : CHR,
#             'SNP' : temp_snp_list1,
#             'BP' : [dic_bp[x] for x in temp_snp_list1],
#             'pCHR' : CHR,
#             'pSNP' : temp_snp_list2,
#             'pBP' : [dic_bp[x] for x in temp_snp_list2],
#             'pAN:proxy_%d_%d_%s_%s' % (LB, UB, term1, term2) : 1
#         })

#         gdreg.util.write_annot(df_pannot, OUT_PATH + "/proxy_%d_%d_%s_%s.chr%d.pannot_hr.gz" % (LB, UB, term1, term2, CHR))
        print("proxy_%d_%d_%s_%s.pannot.gz" % (LB, UB, term1, term2), 'size=%d'% len(temp_snp_list1)) 
    
    print('# Finished, time=%0.1fs'%(time.time() - sys_start_time))
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='gdreg')
    parser.add_argument('--pgen_file', type=str, default=None)
    parser.add_argument('--lb', type=int, default=None)
    parser.add_argument('--ub', type=int, default=None)
    parser.add_argument('--out_path', type=str, default=None)
    
    args = parser.parse_args()
    main(args)