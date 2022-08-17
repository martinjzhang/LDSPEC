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
    
    err_msg = "--pannot needs to be one of ['cS2G_all', 'cS2G_promoter', 'cS2G_other']"
    assert PANNOT in ['cS2G_all', 'cS2G_promoter', 'cS2G_other'], err_msg
    
    # df_snp_chr
    df_snp_chr = gdreg.util.read_pgen(PGEN_FILE)["pvar"]
    df_snp_chr = df_snp_chr[["CHR", "SNP", "BP"]].copy()
    CHR = df_snp_chr["CHR"][0]
    df_snp_chr['MAF'] = gdreg.util.read_pgen(PGEN_FILE)["afreq"]["MAF"].astype(np.float32)
    dic_bp = {x:y for x,y in zip(df_snp_chr["SNP"], df_snp_chr["BP"])}
    print('#########')
    print('df_snp_chr', df_snp_chr.shape)
    print(df_snp_chr.iloc[:3]) 
    
    # df_cs2g
    df_cs2g = pd.read_csv(
        '/n/groups/price/martin/data_GDREG/gene_annotation/cS2G/cS2G_UKBB/cS2G.%s.SGscore.gz' % CHR, sep='\t',
    )
    df_snpmap = pd.read_csv(
        '/n/groups/price/martin/data_GDREG/gene_annotation/cS2G/00_bim/UKBB.%s.info' % CHR, sep=' ',
    )

    temp_dic = {x:y for x,y in zip(df_snpmap['ID'], df_snpmap['RS'])}
    df_cs2g['ID'] = df_cs2g['SNP']
    df_cs2g['SNP'] = [temp_dic[x] for x in df_cs2g['ID']]

    df_cs2g['Link'] = df_cs2g['cS2G'] > 0.5
    df_cs2g['Exon'] = np.array(['Exon' in x for x in df_cs2g['INFO']]) & df_cs2g['Link']
    df_cs2g['Promoter'] = np.array(['Promoter' in x for x in df_cs2g['INFO']]) & df_cs2g['Link']
    df_cs2g['Other'] = (~df_cs2g['Exon']) & (~df_cs2g['Promoter']) & df_cs2g['Link']
    print('#########')
    print('df_cs2g', df_cs2g.shape)
    print('Overlap with df_snp_chr:', len(set(df_snp_chr['SNP']) & set(df_cs2g['SNP'])))
    print(df_cs2g.iloc[:3]) 
    print(df_cs2g[['Link', 'Exon', 'Promoter', 'Other']].sum(axis=0))
    
    # cS2G SNP pairs    
    if PANNOT == 'cS2G_all':
        ind_select = df_cs2g['Link']

    if PANNOT == 'cS2G_promoter':
        ind_select = df_cs2g['Link'] & df_cs2g['Promoter']

    if PANNOT == 'cS2G_other':
        ind_select = df_cs2g['Link'] & df_cs2g['Other']

    temp_df = df_cs2g.loc[ind_select].copy()
    temp_df = temp_df.loc[temp_df['SNP'].isin(df_snp_chr['SNP'])]
    temp_df = temp_df.groupby('GENE').agg({'SNP':list})  

    snp_list1 = []
    snp_list2 = []
    for temp_list in temp_df['SNP']:
        for i in range(len(temp_list)):
            for j in range(i+1, len(temp_list)):
                snp_list1.append(temp_list[i])
                snp_list2.append(temp_list[j])
    
    # Write pannot
    print('pAN:%s' % PANNOT, 'size=%d'% len(snp_list1)) 
    snp_pair_list = [(x,y) for x,y in zip(snp_list1, snp_list2)]
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