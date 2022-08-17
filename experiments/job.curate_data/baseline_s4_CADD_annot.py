import gzip
import pandas as pd
import numpy as np
import time
import argparse

def main(args):
    
    sys_start_time = time.time()
    
    BIM_FILE=args.bim_file
    CADD_FILE=args.cadd_file
    OUT_FILE=args.out_file
    
    df_bim = pd.read_csv(BIM_FILE, sep=' ', header=None)
    df_bim.columns = ["CHR", "SNP", "CM", "BP", "A1", "A2"]
    df_bim.index = ['%d:%d:%s:%s'%(a,b,c,d) for a,b,c,d in zip(df_bim['CHR'],df_bim['BP'],df_bim['A1'],df_bim['A2'])]
    
    df_cadd_raw = pd.read_csv(CADD_FILE, sep='\t', header=1, index_col=None)
    df_cadd_raw.index = [
        '%d:%d:%s:%s'%(a,b,c,d) 
        for a,b,c,d in  zip(df_cadd_raw['#Chrom'],df_cadd_raw['Pos'],df_cadd_raw['Ref'],df_cadd_raw['Alt'])
    ]
    df_cadd_raw.columns = ['CADD:%s'%x for x in df_cadd_raw.columns]
    
    df_cadd = df_bim.copy()
    df_cadd = df_cadd.join(df_cadd_raw)
    df_cadd.to_csv(OUT_FILE, sep="\t", index=False)
    
    print('# Finished, time=%0.1fs'%(time.time() - sys_start_time))
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='gdreg')
    parser.add_argument('--bim_file', type=str, default=None)
    parser.add_argument('--cadd_file', type=str, default=None)
    parser.add_argument('--out_file', type=str, default=None)
    
    args = parser.parse_args()
    main(args)