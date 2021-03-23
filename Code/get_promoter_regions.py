import os
import pandas as pd
import argparse
import multiprocessing as mp
import subprocess
from io import StringIO
import numpy as np
files_path  = "Data/"
parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", help="please provide an input file name containing gene annotations", action = "store")
parser.add_argument("-o","--output", help="please provide an output file name", action = "store")
parser.add_argument("-w","--window",type = int, help="please a bp window around the gene body for defining the promoter region upstream of TSS", action = "store")
parser.add_argument("-m","--multiprocess", help = "select this flag along with all flag to do multiprocessing", action = "store_true")
parser.add_argument("-n", "--cpu",type = int, help = "number of cpu cores to use", action = "store")
args = parser.parse_args()
out_path  = args.output + "_results/"
fl_nm = out_path + args.input + ".csv"
gene_df = pd.read_csv(fl_nm)
gene_df = gene_df.iloc[:,:8]
col_list = list(gene_df)
col_list.remove("transcript_length")
col_list.remove("ensembl_transcript_id")

def get_promoter_region(x):
    gn_df = gene_df.loc[gene_df.ensembl_gene_id == x]
    gn_df = gn_df.sort_values(by =["transcript_length"], ascending = False )
    gn_df1 = pd.DataFrame(gn_df.iloc[0,:]).T
    gn_df2 = gn_df1[col_list]
    strnd = gn_df2.strand.values[0]    
    if strnd == -1:
        end_reg = gn_df2.transcription_start_site.values[0] + int(args.window)
        up_reg = gn_df2.transcription_start_site.values[0]
    else:
        up_reg = gn_df2.transcription_start_site.values[0] - int(args.window)
        end_reg = gn_df2.transcription_start_site.values[0]
    fin_list = [(x,gn_df2.chromosome_name.values[0],up_reg,end_reg)]
    fin_df = pd.DataFrame(fin_list)
    return(fin_df)


out_name = out_path + args.output + "_promoter_regions.csv"
if args.multiprocess:
    pool = mp.Pool(args.cpu)
    list_genes = list(set(gene_df.ensembl_gene_id))
    l_df_l = list(pool.map(get_promoter_region,list_genes))
    l_df = pd.concat(l_df_l)
    l_df.columns = ["ensembl_gene_id", "chromosome_name","promoter_start", "promoter_end"]
    l_df.to_csv(out_name, index = False)
else:
    list_genes = list(set(gene_df.ensembl_gene_id))
    l_df_l == list(map(get_promoter_region,list_genes))
    l_df = pd.concat(l_df_l)
    l_df.columns = ["ensembl_gene_id", "chromosome_name","promoter_start", "promoter_end"]
    l_df.to_csv(out_name, index = False)



