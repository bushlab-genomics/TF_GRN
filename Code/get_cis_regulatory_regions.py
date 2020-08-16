import os
import pandas as pd
import argparse
import multiprocessing as mp
import subprocess
from io import StringIO
import numpy as np
parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", help="please provide an input file name containing gene annotations", action = "store")
parser.add_argument("-o","--output", help="please provide an output file prefix name", action = "store")
parser.add_argument("-w","--window",type = int, help="please a bp window around the gene body for defining the regulatory region(default is 50kB)", action = "store")
parser.add_argument("-c","--ctcf", help="provide the filename for ctcf peaks", action = "store")
parser.add_argument("-m","--multiprocess", help = "select this flag along with all flag to do multiprocessing", action = "store_true")
parser.add_argument("-n", "--cpu",type = int, help = "number of cpu cores to use", action = "store")
args = parser.parse_args()

input_path = "Data/"
out_path = args.output + "_results/"
fl_nm = out_path + args.input + ".csv"
gene_df = pd.read_csv(fl_nm)
gene_df = gene_df.iloc[:,:8]
col_list = list(gene_df)
col_list.remove("transcript_length")
col_list.remove("ensembl_transcript_id")

def get_cis_region_ctcf(x):
    gn_df = gene_df.loc[gene_df.ensembl_gene_id == x]
    gn_df = gn_df.sort_values(by =["transcript_length"], ascending = False )
    gn_df1 = pd.DataFrame(gn_df.iloc[0,:]).T
    gn_df2 = gn_df1[col_list]
    strnd = gn_df2.strand.values[0]
    ctcf_file = input_path + args.ctcf  
    if strnd == -1:
        up_reg = gn_df2.transcription_start_site.values[0] + args.window
        end_reg = gn_df2.start_position.values[0] - args.window
        u_l = (gn_df2.chromosome_name.values[0],gn_df2.transcription_start_site.values[0],  up_reg)
        d_l = (gn_df2.chromosome_name.values[0],end_reg,gn_df2.start_position.values[0])
        u_name = x + "_up.bed"
        d_name = x +"_down.bed"
        udf = pd.DataFrame(u_l).T
        ddf = pd.DataFrame(d_l).T
        udf.to_csv(u_name, sep = "\t", index = False, header = False)
        ddf.to_csv(d_name,sep = "\t", index = False, header = False)
        q1_dwn = "bedtools intersect -b "+  d_name + " -a " + ctcf_file + " -wa -f 1.0"
        q1_up = "bedtools intersect -b " + u_name + " -a " + ctcf_file + " -wa -f 1.0"
        ot_up =  subprocess.check_output(q1_up,shell = True)
        ot1_up = ot_up.decode("utf-8")
        lk_up = StringIO(ot1_up)
        ot_dwn =  subprocess.check_output(q1_dwn,shell = True)
        ot1_dwn = ot_dwn.decode("utf-8")
        lk_dwn = StringIO(ot1_dwn)
        os.system("rm -rf " + u_name)
        os.system("rm -rf " + d_name)
        try:
            up_df = pd.read_csv(lk_up,sep = "\t",header = None)
            mid_u = [up_df[1][i] + np.abs(int((up_df[1][i] - up_df[2][i])/2)) for i in list(up_df.index)]
            max_up_ind = np.argmax([abs(mid_u - gn_df2.transcription_start_site.values[0])])
            max_up = mid_u[max_up_ind]
            up_reg_f_e = max_up
            up_reg_f_s = gn_df2.transcription_start_site.values[0]
        except:
            up_reg_f_s = gn_df2.transcription_start_site.values[0]
            up_reg_f_e = up_reg #args.window
        try:
            down_df = pd.read_csv(lk_dwn,sep = "\t",header = None)
            mid_d = [down_df[2][i] - np.abs(int((down_df[1][i] - down_df[2][i])/2)) for i in list(down_df.index)]
            max_down_ind = np.argmax([abs(mid_d - gn_df2.start_position.values[0])])
            max_down = mid_d[max_down_ind]
            down_reg_f_s = max_down
            down_reg_f_e =  gn_df2.start_position.values[0]
        except:
            down_reg_f_e =  gn_df2.start_position.values[0]
            down_reg_f_s = end_reg 
    else:
        up_reg = gn_df2.transcription_start_site.values[0] - args.window
        end_reg = gn_df2.end_position.values[0] + args.window
        u_l = (gn_df2.chromosome_name.values[0],up_reg,gn_df2.transcription_start_site.values[0])
        d_l = (gn_df2.chromosome_name.values[0],gn_df2.end_position.values[0],end_reg)
        u_name = x + "_up.bed"
        d_name = x +"_down.bed"
        udf = pd.DataFrame(u_l).T
        ddf = pd.DataFrame(d_l).T
        udf.to_csv(u_name, sep = "\t", index = False, header = False)
        ddf.to_csv(d_name,sep = "\t", index = False, header = False)
        q1_dwn = "bedtools intersect -b "+  d_name + " -a " + ctcf_file + " -wa -f 1.0"
        q1_up = "bedtools intersect -b " + u_name + " -a " + ctcf_file + " -wa -f 1.0"
        ot_up =  subprocess.check_output(q1_up,shell = True)
        ot1_up = ot_up.decode("utf-8")
        lk_up = StringIO(ot1_up)
        ot_dwn =  subprocess.check_output(q1_dwn,shell = True)
        ot1_dwn = ot_dwn.decode("utf-8")
        fk_dwn = StringIO(ot1_dwn)
        os.system("rm -rf " + u_name)
        os.system("rm -rf " + d_name)
        try:
            up_df = pd.read_csv(lk_up,sep = "\t",header = None)
            mid_u = [up_df[1][i] + np.abs(int((up_df[1][i] - up_df[2][i])/2)) for i in list(up_df.index)]
            max_up_ind = np.argmax([abs(mid_u - gn_df2.transcription_start_site.values[0])])
            max_up = mid_u[max_up_ind]
            up_reg_f_s = max_up
            up_reg_f_e = gn_df2.transcription_start_site.values[0]
        except:
            up_reg_f_e = gn_df2.transcription_start_site.values[0]
            up_reg_f_s = up_reg
        try:    
            down_df = pd.read_csv(lk_dwn,sep = "\t",header = None)
            mid_d = [down_df[2][i] - np.abs(int((down_df[1][i] - down_df[2][i])/2)) for i in list(down_df.index)]
            max_down_ind = np.argmax([abs(mid_d - gn_df2.end_position.values[0])])
            max_down = mid_d[max_down_ind]
            down_reg_f_e = max_down
            down_reg_f_s =  gn_df2.end_position.values[0]
        except:
           down_reg_f_s =  gn_df2.end_position.values[0]
           down_reg_f_e = end_reg
    fin_list = [(x,gn_df2.chromosome_name.values[0],up_reg_f_s,up_reg_f_e),(x,gn_df2.chromosome_name.values[0],down_reg_f_s,down_reg_f_e)]
    fin_df = pd.DataFrame(fin_list)
    return(fin_df)


def get_cis_region(x):
    gn_df = gene_df.loc[gene_df.ensembl_gene_id == x]         
    gn_df = gn_df.sort_values(by =["transcript_length"], ascending = False )
    gn_df1 = pd.DataFrame(gn_df.iloc[0,:]).T
    gn_df2 = gn_df1[col_list]
    strnd = gn_df2.strand.values[0]
    if strnd == -1:
        up_reg = gn_df2.transcription_start_site.values[0] + args.window
        end_reg = gn_df2.start_position.values[0] - args.window
        down_reg_f_e =  gn_df2.start_position.values[0]
        down_reg_f_s = end_reg
        up_reg_f_s = gn_df2.transcription_start_site.values[0]
        up_reg_f_e = up_reg #args.window
    else:
        up_reg = gn_df2.transcription_start_site.values[0] - args.window
        end_reg = gn_df2.end_position.values[0] + args.window
        down_reg_f_s =  gn_df2.end_position.values[0]    
        down_reg_f_e = end_reg 
        up_reg_f_e = gn_df2.transcription_start_site.values[0]
        up_reg_f_s = up_reg 
    fin_list = [(x,gn_df2.chromosome_name.values[0],up_reg_f_s,up_reg_f_e),(x,gn_df2.chromosome_name.values[0],down_reg_f_s,down_reg_f_e)]
    fin_df = pd.DataFrame(fin_list)
    return(fin_df)


def robust_cis_ctcf(y):
        try:
                dft = get_cis_region_ctcf(y)
        except:
                dct = [y,np.nan,np.nan,np.nan]
                dft = pd.DataFrame(dct).T
        return(dft)

def robust_cis(y):
        try:
                dft = get_cis_region(y)        
        except:
                dct = [y,np.nan,np.nan,np.nan]                                                                                                
                dft = pd.DataFrame(dct).T
        return(dft)


out_name = out_path +  args.output + "_reg_regions.csv"
if args.ctcf != "none":
    if args.multiprocess:
        pool = mp.Pool(args.cpu)
        list_genes = sorted(set(gene_df.ensembl_gene_id))
        l_df_l = list(pool.map(robust_cis_ctcf,list_genes))
        #l_df_l = list(pool.map(get_cis_region_ctcf,list_genes))
        l_df = pd.concat(l_df_l)
        l_df.columns = ["ensembl_gene_id", "chromosome_name", "reg_start", "reg_end"]
        l_df = l_df.dropna()
        l_df.to_csv(out_name, index = False)
    else:
        list_genes = sorted(set(gene_df.ensembl_gene_id))
        l_df_l = list(map(robust_cis_ctcf,list_genes))       
        l_df = pd.concat(l_df_l)
        l_df.columns = ["ensembl_gene_id", "chromosome_name", "reg_start", "reg_end"]
        l_df = l_df.dropna()
        l_df.to_csv(out_name, index = False)
else:
    if args.multiprocess:
        pool = mp.Pool(args.cpu)
        list_genes = list(set(gene_df.ensembl_gene_id))
        l_df_l = list(pool.map(robust_cis,list_genes))
        l_df = pd.concat(l_df_l)
        l_df.columns = ["ensembl_gene_id", "chromosome_name", "reg_start", "reg_end"]
        l_df = l_df.dropna()
        l_df.to_csv(out_name, index = False)
    else:
        list_genes = list(set(gene_df.ensembl_gene_id))
        l_df_l = list(map(robust_cis,list_genes))
        l_df = pd.concat(l_df_l)
        l_df.columns = ["ensembl_gene_id", "chromosome_name", "reg_start", "reg_end"]
        l_df = l_df.dropna()
        l_df.to_csv(out_name, index = False)

os.system("rm -rf ENS*.bed")


