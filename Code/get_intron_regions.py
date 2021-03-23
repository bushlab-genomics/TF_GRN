import pandas as pd
import os
import multiprocessing as mp
import subprocess
from io import StringIO
import glob as glob
#files_path = "Data/"
import argparse


parser = argparse.ArgumentParser()

parser.add_argument("-i","--input", help="please provide a file name containing input file name containing gene exon info", action = "store")
parser.add_argument("-o","--output", help="please provide an output file name", action = "store")
parser.add_argument("-r", "--reg", help= "please provide a file name containing gene regulatory regions", action = "store")
parser.add_argument("-m","--multiprocess", help = "select this flag along with all flag to do multiprocessing", action = "store_true")
parser.add_argument("-n", "--cpu",type = int, help = "number of cpu cores to use", action = "store")

args = parser.parse_args()
out_path  = args.output + "_results/"
t_e_df = pd.read_csv(out_path + args.input)
reg_df = pd.read_csv(out_path + args.reg)

def run_t_bedtools(y):
	tr_df = t_e_df.loc[t_e_df.ensembl_gene_id == y]
	tr_df1 = tr_df[["chromosome_name","transcript_start","transcript_end"]]
	tr_df1 = tr_df1.drop_duplicates()
	tr_df2 =  tr_df[["chromosome_name", "exon_chrom_start", "exon_chrom_end"]]
	tr_df_nm = out_path + y + "_transcript.bed"
	e_df_nm = out_path +y +"_exon.bed"
	tr_df1.to_csv(tr_df_nm, sep = "\t", header = False, index = False)
	tr_df2.to_csv(e_df_nm, sep = "\t", header = False, index = False)
	g_bed_nm = out_path + y + ".bed"			
	q1 = "bedtools subtract -a " + tr_df_nm  + " -b " + e_df_nm + " > " +  g_bed_nm 
	os.system(q1)
	kl = pd.read_csv(g_bed_nm,sep = "\t",header = None)
	q2 = "rm " + tr_df_nm
	q3 = "rm " + e_df_nm
	q4 = "rm " + g_bed_nm
	os.system(q2)
	os.system(q3)
	os.system(q4)
	kl.columns  = ["chromosome_name", "reg_start", "reg_end"]
	kl.insert(0, "ensembl_gene_id", [y]*len(kl))
	return(kl)
	


def robust_run(y):
	try:
		dft = run_t_bedtools(y)
	except:
		dct = {"ensembl_gene_id":[y], "chromosome_name":["NA"],"reg_start":["NA"], "reg_end":["NA"]}
		dft = pd.DataFrame(dct)
	return(dft)


if args.multiprocessing:
	pool = mp.Pool(args.cpu)
	gene_list = list(set(t_e_df.ensembl_gene_id))
	tdf_list = list(pool.map(robust_run, gene_list))
else:
	gene_list = list(set(t_e_df.ensembl_gene_id))
	tdf_list = list(map(robust_run, gene_list))
tdf = pd.concat(tdf_list)
tdf = tdf.dropna()
df_f = pd.concat([tdf, reg_df])
df_f = df_f.dropna()
o_name =  out_path + args.output + "_intron_regulatory_regions.csv"
df_f.to_csv(o_name, index = False)
os.system("rm " + out_path + "ENS*.bed")



