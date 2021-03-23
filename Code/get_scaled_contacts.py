import pandas as pd
import os
import pandas as pd
import argparse
#import multiprocessing as mp
import subprocess
from io import StringIO
import numpy as np
from sklearn.preprocessing import MinMaxScaler,RobustScaler

#hic_path  = "Data/"
parser = argparse.ArgumentParser()
parser.add_argument("-t","--tfbs", help="please provide a filename for all tfbs", action = "store")
parser.add_argument("-c","--hic", help="please provide a filename for hic contacts", action = "store")
parser.add_argument("-p","--promoter", help="please provide a filename for promoter tfbs", action = "store") 
parser.add_argument("-s","--scalingfactor",  help="please provide a scaling factor in the range [0,1]", action = "store")
parser.add_argument("-u","--upweighting",  help="select this flag if you want  to upweight the promoter TFBS in the HiC analysis", action = "store_true")
parser.add_argument("-o","--output",  help="please provide an output file name", action = "store")

args = parser.parse_args() 
path = args.output + "_results/" 
all_tfbs = pd.read_csv(path + args.tfbs) 
prom_df = pd.read_csv(path + args.promoter) 
hic_df = pd.read_csv(path + args.hic) 
col_list = ["Ensembl","Chromosome", "TF_Start", "TF_End","TF_Name"] 
hic_df1 = hic_df[["Ensembl","Chromosome", "TF_Start", "TF_End","TF_Name", "Contacts"]] 
sf = args.scalingfactor 
np_contacts = np.array(hic_df1.Contacts) 
scaled_contacts = MinMaxScaler().fit_transform(np_contacts.reshape((-1,1))) + float(sf) 
hic_df1["Contacts"] = scaled_contacts 
prom_df["prom_index"] = list(prom_df.index) 
all_tfbs["reg_index"] = list(all_tfbs.index) 
hic_prom_merge = hic_df1.merge(prom_df, on =col_list, how = "inner") 
all_tfbs_hic_merge = all_tfbs.merge(hic_df1, on = col_list , how = "inner") 
tfbs_promoter_merge = all_tfbs.merge(prom_df, on = col_list , how = "inner") 
if args.upweighting:
	sf = 2.0
else:
	sf = 1.0

if len(hic_prom_merge) == 0:
	df1 = hic_df1
	df2 = all_tfbs.loc[~all_tfbs.reg_index.isin(all_tfbs_hic_merge.reg_index)]
	df2 = df2.loc[~df2.reg_index.isin(tfbs_promoter_merge.reg_index)]
	df2 = df2[col_list]
	df2["Contacts"] = [1.0]*len(df2)
	df3 = all_tfbs.loc[all_tfbs.reg_index.isin(tfbs_promoter_merge.reg_index)]
	df3 = df3[col_list]
	df3["Contacts"] = [sf]*len(df3)
	df_f = pd.concat([df1,df2,df3])
else:
	df1 = hic_prom_merge
	df1 = df1[["Ensembl","Chromosome", "TF_Start", "TF_End","TF_Name", "Contacts"]]
	df2  =  all_tfbs.loc[~all_tfbs.reg_index.isin(all_tfbs_hic_merge.reg_index)]
	df2 = df2.loc[~df2.reg_index.isin(tfbs_promoter_merge.reg_index)]
	df2 = df2[col_list]
	df2["Contacts"] = [1.0]*len(df2)
	df3 = prom_df.loc[~prom_df.prom_index.isin(hic_prom_merge.prom_index)]
	df4 = all_tfbs.merge(df3, on = col_list, how = "inner")
	df4 = df4[col_list]
	df4["Contacts"] = [sf]*len(df4)	
	df_f = pd.concat([df1, df2, df4])

df_f0 = df_f[["Ensembl", "TF_Name", "Contacts"]]
df_f1 = df_f0.groupby(["TF_Name", "Ensembl"]).mean()
ind_list = [[i[0],i[1]] for i in list(df_f1.index)]
ind_df = pd.DataFrame(ind_list)
ind_df.columns= ["TF_Name", "Ensembl"]
df_f1 = df_f1.reset_index(drop = True)
df_f2 = pd.concat([ind_df, df_f1], axis=1)
o_name =  path + args.output + "_hic_tfbs.csv"
df_f2.to_csv(o_name, index = False)





