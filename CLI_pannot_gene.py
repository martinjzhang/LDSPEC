import gzip
import pandas as pd
import numpy as np
import time
import argparse
import ldspec

"""
Create gene-based pannot annotations.

Basic gene-based annotation: 
    --pgen_file 
    --pannot
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
    PANNOT = args.pannot
    LD_FILE = args.ld_file
    SNP_RANGE_FILE = args.snp_range_file
    LD_LB = float(args.ld_lb)
    LD_UB = float(args.ld_ub)
    MAF_BIN_FILE = args.maf_bin_file
    MAF_RATIO_THRES = int(args.maf_ratio_thres)
    FLAG_MAF_BLOCK = args.flag_maf_block
    OUT_PATH = args.out_path

    print("--pgen_file %s" % PGEN_FILE)
    print("--pannot %s" % PANNOT)
    print("--ld_file %s" % LD_FILE)
    print("--snp_range_file %s" % SNP_RANGE_FILE)
    print("--ld_lb %s" % LD_LB)
    print("--ld_ub %s" % LD_UB)
    print("--maf_bin_file %s" % MAF_BIN_FILE)
    print("--maf_ratio_thres %s" % MAF_RATIO_THRES)
    print("--flag_maf_block %s" % FLAG_MAF_BLOCK)
    print("--out_path %s" % OUT_PATH)
    
    err_msg = "--pannot needs to be one of ['gene', 'exon', 'exonic_gene', 'protein_domain', 'cS2G_promoter']"
    assert PANNOT in ['gene', 'exon', 'exonic_gene', 'protein_domain', 'cS2G_promoter'], err_msg

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
    
    # Data loading: gene info
    if PANNOT in ['gene', 'exon', 'exonic_gene', 'protein_domain']:
        df_gene = pd.read_csv("/n/groups/price/martin/LDSPEC_data/gene_annotation/ENSG_gene_annot_v41.txt", sep="\t")
        df_exon = pd.read_csv("/n/groups/price/martin/LDSPEC_data/gene_annotation/ENSE_exon_annot_v41.txt", sep="\t")
        df_gene_chr = df_gene.loc[df_gene["CHR"] == 'chr%d'%CHR]
        df_exon_chr = df_exon.loc[df_exon["CHR"] == 'chr%d'%CHR]     
        v_bp = df_snp_chr['BP'].values
        print(df_gene_chr.iloc[:3])
        print(df_exon_chr.iloc[:3])
    if PANNOT in ['cS2G_promoter']:
        df_cs2g = pd.read_csv(
        '/n/groups/price/martin/LDSPEC_data/gene_annotation/cS2G/cS2G_UKBB/cS2G.%s.SGscore.gz' % CHR, sep='\t',
        )
        df_snpmap = pd.read_csv(
            '/n/groups/price/martin/LDSPEC_data/gene_annotation/cS2G/00_bim/UKBB.%s.info' % CHR, sep=' ',
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

    # Computation
    snp_list = list(df_snp_chr["SNP"])
    snp_list1 = []
    snp_list2 = []
    dic_pair_gene = {x:set() for x in snp_list} # dic_pair_gene[snp] = (snp1, snp2, ...)

    if PANNOT == 'gene':
        for START,END in zip(df_gene_chr["START"], df_gene_chr["END"]):
            ind_select = (v_bp>=START) & (v_bp<=END)
            temp_list = df_snp_chr["SNP"].values[ind_select]
            for i in range(len(temp_list)):
                for j in range(i+1, len(temp_list)):
                    dic_pair_gene[temp_list[i]].add(temp_list[j])
    if PANNOT == 'exon':
        for START,END in zip(df_exon_chr["START"], df_exon_chr["END"]):
            ind_select = (v_bp>=START) & (v_bp<=END)
            temp_list = df_snp_chr["SNP"].values[ind_select]
            for i in range(len(temp_list)):
                for j in range(i+1, len(temp_list)):
                    dic_pair_gene[temp_list[i]].add(temp_list[j])
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
                    dic_pair_gene[temp_list[i]].add(temp_list[j])                   
    if PANNOT == 'protein_domain':
        file_folder = '/n/groups/price/martin/LDSPEC_data/UKBimp_337K_MAF001/baseline_annot/vep'
        df_pd = ldspec.util.read_annot(file_folder + '/ukb_imp_chr%s_v3.vep.annot.gz' % CHR)
        dic_pd = {x:set(y.split(',')) for x,y in zip(df_pd['SNP'], df_pd['AN:DOMAINS'])}

        v_flag_pd = np.array([x in dic_pd for x in df_snp_chr["SNP"]], dtype=bool)
        for START,END in zip(df_gene_chr["START"], df_gene_chr["END"]):
            ind_select = (v_bp>=START) & (v_bp<=END) & v_flag_pd
            temp_list = df_snp_chr["SNP"].values[ind_select]
            for i in range(len(temp_list)):
                for j in range(i+1, len(temp_list)):
                    if len(dic_pd[temp_list[i]] & dic_pd[temp_list[j]]) > 0:
                        dic_pair_gene[temp_list[i]].add(temp_list[j])
    if PANNOT == 'cS2G_promoter':
        ind_select = df_cs2g['Link'] & df_cs2g['Promoter']
        temp_df = df_cs2g.loc[ind_select].copy()
        temp_df = temp_df.loc[temp_df['SNP'].isin(df_snp_chr['SNP'])]
        temp_df = temp_df.groupby('GENE').agg({'SNP':list})          
        for temp_list in temp_df['SNP']:
            for i in range(len(temp_list)):
                for j in range(i+1, len(temp_list)):
                    dic_pair_gene[temp_list[i]].add(temp_list[j])
                    
    # Filter by LD
    if LD_FILE is None:       
        for snp1 in dic_pair_gene:    
            for snp2 in dic_pair_gene[snp1]:
                snp_list1.append(snp1)
                snp_list2.append(snp2)
    else:
        for snp_range in snp_range_list_chr:
            print("    snp_range=%s" % snp_range)
            mat_ld, dic_range = ldspec.util.read_ld(LD_FILE.replace("@", snp_range))        
            for i in range(dic_range["start"], dic_range["end"]):        
                if snp_list[i] not in dic_pair_gene:
                    continue
                v_ld = mat_ld[:, i - dic_range["start"]].toarray().flatten()
                temp_set = dic_pair_gene[snp_list[i]]
                temp_list = [x for x in range(i+1, len(snp_list)) if snp_list[x] in temp_set]
                for j in temp_list:
                    ld = v_ld[j]
                    if (ld >= LD_LB) & (ld <= LD_UB):
                        snp_list1.append(snp_list[i])
                        snp_list2.append(snp_list[j])

    print("n_pair", len(snp_list1))

    # MAF-bin & write files
    LD_LB = max(-1, LD_LB)
    LD_UB = min(1, LD_UB)
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
            file_name = "%s.%s.%s.c%d" % (PANNOT, str_ld, str_maf, CHR)
            print("%-50s" % file_name, "n_pair=%d" % len(snp_pair_list))
            ldspec.util.write_pannot_mat(
                snp_pair_list, list(df_snp_chr["SNP"]), OUT_PATH + "/" + file_name
            )

    print("# Finished, time=%0.1fs" % (time.time() - sys_start_time))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ldspec")
    parser.add_argument("--pgen_file", type=str, default=None, required=True)
    parser.add_argument('--pannot', type=str, default=None, required=True)    
    parser.add_argument("--ld_file", type=str, default=None, required=False)
    parser.add_argument("--snp_range_file", type=str, default=None, required=False)
    parser.add_argument("--ld_lb", type=str, default="-1", required=False)
    parser.add_argument("--ld_ub", type=str, default="1", required=False)
    parser.add_argument("--maf_bin_file", type=str, default=None, required=False)
    parser.add_argument("--maf_ratio_thres", type=str, default="5", required=False)
    parser.add_argument("--flag_maf_block", type=bool, default=True, required=False)
    parser.add_argument("--out_path", type=str, default=None, required=True)

    args = parser.parse_args()
    main(args)
