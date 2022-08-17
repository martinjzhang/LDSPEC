import gzip
import pandas as pd
import numpy as np
import time
import argparse
import gdreg

def main(args):
    
    sys_start_time = time.time()
    
    PGEN_FILE = args.pgen_file
    OUT_PATH = args.out_path
    
    print("exon .pannot_hr.gz")
    print("--pgen_file %s" % PGEN_FILE)
    print("--out_path %s" % OUT_PATH)
    
    # Exon .pannot_mat.gz
    GENE_FILE = "/n/groups/price/martin/data_GDREG/gene_annotation/ENSE_exon_annot_v41.txt"
    df_gene = pd.read_csv(GENE_FILE, sep="\t")
    
    df_snp_chr = gdreg.util.read_pgen(PGEN_FILE)["pvar"]
    df_snp_chr = df_snp_chr[["CHR", "SNP", "BP"]].copy()
    CHR = df_snp_chr["CHR"][0]
    df_snp_chr['MAF'] = gdreg.util.read_pgen(PGEN_FILE)["afreq"]["MAF"].astype(np.float32)
    dic_bp = {x:y for x,y in zip(df_snp_chr["SNP"], df_snp_chr["BP"])}

    # Exon-level SNP pairs
    temp_df = df_gene.loc[df_gene["CHR"] == 'chr%d'%CHR]
    snp_list1 = []
    snp_list2 = []
    for START,END in zip(temp_df["START"], temp_df["END"]):
        temp_list = df_snp_chr["SNP"].values[(df_snp_chr['BP']>=START) & (df_snp_chr['BP']<=END)]
        for i in range(len(temp_list)):
            for j in range(i+1, len(temp_list)):
                snp_list1.append(temp_list[i])
                snp_list2.append(temp_list[j])

    dic_maf_snp = {
        "common" : set(df_snp_chr["SNP"][df_snp_chr["MAF"].values >= 0.05]),
        "lf" : set(df_snp_chr["SNP"][(df_snp_chr["MAF"].values >= 0.005) & (df_snp_chr["MAF"].values < 0.05)]),
        "rare" : set(df_snp_chr["SNP"][df_snp_chr["MAF"].values <= 0.005]),
    }

    for term1,term2 in [
        ['common', 'common'], 
        ['common', 'lf'], 
        ['common', 'rare'], 
        ['lf', 'lf'], 
        ['lf', 'rare'],
        ['rare', 'rare']
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
                snp_pair_list, list(df_snp_chr["SNP"]), OUT_PATH + "/exon_%s_%s.chr%d" % (term1, term2, CHR)
            )
        print('pAN:exon_%s_%s' % (term1, term2), 'size=%d'% len(temp_snp_list1)) 
    
    print('# Finished, time=%0.1fs'%(time.time() - sys_start_time))
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='gdreg')
    parser.add_argument('--pgen_file', type=str, default=None)
    parser.add_argument('--out_path', type=str, default=None)
    
    args = parser.parse_args()
    main(args)