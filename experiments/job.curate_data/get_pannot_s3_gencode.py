import gzip
import pandas as pd
import numpy as np
import time
import argparse
import gdreg

def main(args):
    
    sys_start_time = time.time()
    
    PGEN_FILE = args.pgen_file
    PANNOT = args.pannot
    OUT_PATH = args.out_path
    
    print("--pgen_file %s" % PGEN_FILE)
    print("--pannot %s" % PANNOT)
    print("--out_path %s" % OUT_PATH)
    
    err_msg = "--pannot needs to be one of ['gene', 'exon', 'exonic_gene', 'protein_domain']"
    assert PANNOT in ['gene', 'exon', 'exonic_gene', 'protein_domain'], err_msg
    
    df_gene = pd.read_csv("/n/groups/price/martin/data_GDREG/gene_annotation/ENSG_gene_annot_v41.txt", sep="\t")
    print('#########')
    print('df_gene', df_gene.shape)
    print(df_gene.iloc[:3])
    df_exon = pd.read_csv("/n/groups/price/martin/data_GDREG/gene_annotation/ENSE_exon_annot_v41.txt", sep="\t")
    print('#########')
    print('df_exon', df_exon.shape)
    print(df_exon.iloc[:3])
    
    # df_snp_chr
    df_snp_chr = gdreg.util.read_pgen(PGEN_FILE)["pvar"]
    df_snp_chr = df_snp_chr[["CHR", "SNP", "BP"]].copy()
    CHR = df_snp_chr["CHR"][0]
    df_snp_chr['MAF'] = gdreg.util.read_pgen(PGEN_FILE)["afreq"]["MAF"].astype(np.float32)
    print('#########')
    print('df_snp_chr', df_snp_chr.shape)
    print(df_snp_chr.iloc[:3]) 

    # Gene-level SNP pairs
    df_gene_chr = df_gene.loc[df_gene["CHR"] == 'chr%d'%CHR]
    df_exon_chr = df_exon.loc[df_exon["CHR"] == 'chr%d'%CHR]    
    v_bp = df_snp_chr['BP'].values       
    snp_list1 = []
    snp_list2 = []            
    
    if PANNOT == 'gene':
        for START,END in zip(df_gene_chr["START"], df_gene_chr["END"]):
            ind_select = (v_bp>=START) & (v_bp<=END)
            temp_list = df_snp_chr["SNP"].values[ind_select]
            for i in range(len(temp_list)):
                for j in range(i+1, len(temp_list)):
                    snp_list1.append(temp_list[i])
                    snp_list2.append(temp_list[j])
                    
    if PANNOT == 'exon':
        for START,END in zip(df_exon_chr["START"], df_exon_chr["END"]):
            ind_select = (v_bp>=START) & (v_bp<=END)
            temp_list = df_snp_chr["SNP"].values[ind_select]
            for i in range(len(temp_list)):
                for j in range(i+1, len(temp_list)):
                    snp_list1.append(temp_list[i])
                    snp_list2.append(temp_list[j])
                    
    if PANNOT == 'exonic_gene':
        # Flag for if the SNP is in exon
        v_flag_exon = np.zeros(df_snp_chr.shape[0], dtype=bool)
        for START,END in zip(df_exon_chr["START"], df_exon_chr["END"]):
            ind_select = (v_bp>=START) & (v_bp<=END)
            v_flag_exon[ind_select] = True
            
        for START,END in zip(df_gene_chr["START"], df_gene_chr["END"]):
            ind_select = (v_bp>=START) & (v_bp<=END) & v_flag_exon
            temp_list = df_snp_chr["SNP"].values[ind_select]
            for i in range(len(temp_list)):
                for j in range(i+1, len(temp_list)):
                    snp_list1.append(temp_list[i])
                    snp_list2.append(temp_list[j])
                    
    if PANNOT == 'protein_domain':
        file_folder = '/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/baseline_annot/vep'
        df_pd = gdreg.util.read_annot(file_folder + '/ukb_imp_chr%s_v3.vep.annot.gz' % CHR)
        dic_pd = {x:set(y.split(',')) for x,y in zip(df_pd['SNP'], df_pd['AN:DOMAINS'])}

        v_flag_pd = np.array([x in dic_pd for x in df_snp_chr["SNP"]], dtype=bool)
        for START,END in zip(df_gene_chr["START"], df_gene_chr["END"]):
            ind_select = (v_bp>=START) & (v_bp<=END) & v_flag_pd
            temp_list = df_snp_chr["SNP"].values[ind_select]
            for i in range(len(temp_list)):
                for j in range(i+1, len(temp_list)):
                    if len(dic_pd[temp_list[i]] & dic_pd[temp_list[j]]) > 0:
                        snp_list1.append(temp_list[i])
                        snp_list2.append(temp_list[j])
    
    # Write pannot
    print('pAN:%s' % PANNOT, 'size=%d'% len(snp_list1)) 
    snp_pair_list = [(x,y) for x,y in zip(snp_list1, snp_list1)]
    if len(snp_pair_list) > 10:
        gdreg.util.write_pannot_mat(
            snp_pair_list, list(df_snp_chr["SNP"]), OUT_PATH + "/%s.chr%d" % (PANNOT, CHR)
        )

    # Stratify by MAF bins and 
    dic_maf_snp = {
        "common" : set(df_snp_chr["SNP"][df_snp_chr["MAF"].values >= 0.05]),
        "lf" : set(df_snp_chr["SNP"][(df_snp_chr["MAF"].values >= 0.005) & (df_snp_chr["MAF"].values < 0.05)]),
        "rare" : set(df_snp_chr["SNP"][df_snp_chr["MAF"].values <= 0.005]),
    }
    
    for mbin1,mbin2 in [['common', 'common'], ['common', 'lf'], ['lf', 'lf']]:
        temp_snp_list1 = []
        temp_snp_list2 = []
        for snp1,snp2 in zip(snp_list1, snp_list2):
            flag1 = (snp1 in dic_maf_snp[mbin1]) & (snp2 in dic_maf_snp[mbin2])
            flag2 = (snp1 in dic_maf_snp[mbin2]) & (snp2 in dic_maf_snp[mbin1])
            if flag1 | flag2:
                temp_snp_list1.append(snp1)
                temp_snp_list2.append(snp2)                  
                
        snp_pair_list = [(x,y) for x,y in zip(temp_snp_list1, temp_snp_list2)]
        if len(snp_pair_list) > 10:
            gdreg.util.write_pannot_mat(
                snp_pair_list, list(df_snp_chr["SNP"]), OUT_PATH + "/%s_%s_%s.chr%d" % (PANNOT, mbin1, mbin2, CHR)
            )
        print('pAN:%s_%s_%s' % (PANNOT, mbin1, mbin2), 'size=%d'% len(temp_snp_list1)) 
    
    print('# Finished, time=%0.1fs'%(time.time() - sys_start_time))
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='gdreg')
    parser.add_argument('--pgen_file', type=str, default=None)
    parser.add_argument('--pannot', type=str, default=None)    
    parser.add_argument('--out_path', type=str, default=None)
    
    args = parser.parse_args()
    main(args)