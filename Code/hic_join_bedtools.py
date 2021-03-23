import numpy as np
import pandas as pd
import os
import multiprocessing as mp
import subprocess
from io import StringIO
import itertools as itr
import argparse 

parser = argparse.ArgumentParser()
parser.add_argument("-p","--promoter", help="please provide a file containing promter regions", action = "store")
parser.add_argument("-t","--tfbs", help="please provide a file containing all tfbs", action = "store")
parser.add_argument("-c","--hic", help="please provide a file containing hic contacts", action = "store")
parser.add_argument("-o","--output", help="please provide a filename for the output file", action = "store")
args = parser.parse_args()

path = args.output + "_results/"
hic_path  =  "Data/"
prom_df = pd.read_csv(path + args.promoter)
tfbs_df = pd.read_csv(path + args.tfbs)
hic_df = pd.read_csv(hic_path + args.hic)
prom_df.columns =  ["Ensembl", "Chromosome", "promoter_start", "promoter_end"]
prom_df["Chromosome"]= prom_df["Chromosome"].astype("int")
tp_merge = prom_df.merge(tfbs_df,  on =["Ensembl", "Chromosome"], how = "inner")
tp_merge = tp_merge.drop_duplicates()
tp_merge = tp_merge.dropna()
tp_merge["reg_index"] = list(tp_merge.index)

def parse_sort_df(df):
    df1 = df.iloc[:,[0,1,6,7]]
    df2 = df.iloc[:,[2,3,4,5]]
    s_list = [sorted([sorted((list(df2.iloc[i,[0,1]]))), sorted(list(df2.iloc[i,[2,3]]))]) for i in list(df2.index)]
    s_list1 = [list(itr.chain.from_iterable(i)) for i in s_list]
    s_df = pd.DataFrame(s_list1)
    s_df.columns = ["up_start", "up_end", "down_start", "down_end"]
    df3 = pd.concat([df1,s_df], axis = 1)
    #df3["reg_index"] = list(df3.index)
    return(df3)


tp_sort_df = parse_sort_df(tp_merge)
tp_up_reg = tp_sort_df.iloc[:,[1,4,5,3]]
tp_down_reg = tp_sort_df.iloc[:,[1,6,7,3]]
tp_up_reg.to_csv(path + "up_reg_temp.bed", sep = "\t", header = None, index = False)
tp_down_reg.to_csv(path + "down_reg_temp.bed", sep = "\t", header =None, index = False)
hic_df["cont_index"] = list(hic_df.index)

up_hic_df = hic_df[["Chromosome", "Start","Start","Contacts", "cont_index"]]
down_hic_df = hic_df[["Chromosome", "Stop","Stop", "Contacts","cont_index"]]

up_hic_df.to_csv(path + "up_cont_temp.bed", sep = "\t", header = None, index = False)
down_hic_df.to_csv(path + "down_cont_temp.bed",sep = "\t", header =None, index = False)

q_up = "bedtools intersect -a " + path + "up_reg_temp.bed -b " + path + "up_cont_temp.bed -wa -wb > " + path + "up_merge_temp.bed"
q_down = "bedtools intersect -a " +  path  + "down_reg_temp.bed -b " + path  + "down_cont_temp.bed -wa -wb > " + path +  "down_merge_temp.bed"
os.system(q_up)
os.system(q_down)

u_merge_df = pd.read_csv(path + "up_merge_temp.bed", sep ="\t", header = None)
d_merge_df = pd.read_csv(path + "down_merge_temp.bed", sep = "\t", header = None)

u_merge_df = u_merge_df.iloc[:,[0,1,2,3,5,7,8]]
d_merge_df = d_merge_df.iloc[:,[0,1,2,3,5,8]]

u_merge_df.columns = ["Chromosome", "up_start","up_end", "reg_index", "Start", "Contacts", "cont_index"]
d_merge_df.columns = ["Chromosome", "down_start","down_end", "reg_index", "Stop","cont_index"]

all_merge_df = u_merge_df.merge(d_merge_df, on = ["Chromosome", "reg_index", "cont_index"], how = "inner")
all_merge_df = all_merge_df[["reg_index", "Start", "Stop","Contacts"]]

#all_merge_df1 = tp_sort_df.merge(all_merge_df, on = ["Chromosome", "up_start", "up_end", "down_start", "down_end", "reg_index"], how = "inner")
all_merge_df1 = tp_merge.merge(all_merge_df, on = ["reg_index"], how = "inner") 
o_name = path + args.output + "_hic_merge.csv"
all_merge_df1.to_csv(o_name, index = False)
q_rm  = "rm " + path + "*temp*"
os.system(q_rm)








