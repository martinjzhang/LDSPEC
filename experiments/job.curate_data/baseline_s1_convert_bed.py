import pandas as pd
import numpy as np
import os
import time
import argparse
import gdreg

def main(args):
    
    sys_start_time = time.time()
    BED_FILE_PATH="/n/groups/price/ldsc/reference_files/bed"
    
    PGEN_FILE=args.pgen_file
    PREFIX_OUT=args.prefix_out    
    print("PGEN_FILE=%s" % PGEN_FILE)
    print("PREFIX_OUT=%s" % PREFIX_OUT)
    
    BED_FILE_LIST = []
    temp_list = os.listdir(BED_FILE_PATH)
    with open("bed_list.txt", "r") as f:
        for line in f:
            line = line.strip()
            if "%s.bed" % line in temp_list:
                BED_FILE_LIST.append("%s.bed" % line)
            if "%s.extend.500.bed" % line in temp_list:
                BED_FILE_LIST.append("%s.extend.500.bed" % line)
    BED_FILE_LIST = sorted(BED_FILE_LIST)
    print("Detected %d .bed files " % len(BED_FILE_LIST))
    
    # Read PGEN
    dic_data = gdreg.util.read_pgen(PGEN_FILE)
    CHR = dic_data["pvar"]['CHR'][0]
    
    # Create df_annot 
    df_annot = dic_data["pvar"][["CHR", "SNP", "BP"]].copy()

    for bed_file in BED_FILE_LIST:
        
        fpath = BED_FILE_PATH + "/" + bed_file
        df_bed = pd.read_csv(fpath, sep='\t', header=None)
        df_bed = df_bed.loc[df_bed[0]=='chr%d'%CHR]
        if df_bed.shape[1] == 3:
            df_bed[3] = 1
        df_bed.columns = ['CHR', "START", "END", "VAL"]
        df_bed["VAL"] = df_bed["VAL"].astype(np.float32)

        dic_bed = {}
        for START,END,VAL in zip(df_bed["START"], df_bed["END"], df_bed["VAL"]):
            # pvar files seem to be 1-based 
            # bed files are 0-based
            temp_dic = {x:VAL for x in range(START+1, END+1)}
            dic_bed.update(temp_dic)

        bed_name = "AN:" + bed_file.replace('.bed', '')
        df_annot[bed_name] = [dic_bed[x] if x in dic_bed else 0 for x in df_annot['BP']]

        if len(df_annot[bed_name].unique())==2:
            df_annot[bed_name] = df_annot[bed_name].astype(bool)
        else:
            df_annot[bed_name] = df_annot[bed_name].astype(np.float32)

        print("%s : %d unique values, %s non-zeros" % (
            bed_name, len(df_annot[bed_name].unique()), (df_annot[bed_name]!=0).sum()
        ))
    
    # Write df_annot
    gdreg.util.write_annot(df_annot, PREFIX_OUT+'.bed.annot.gz')
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='gdreg')
    parser.add_argument('--pgen_file', type=str, default=None)
    parser.add_argument('--prefix_out', type=str, default=None)
    
    args = parser.parse_args()
    main(args)