import os
import pandas as pd
import argparse
import multiprocessing as mp
import subprocess
from io import StringIO
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-f","--file", help="please provide a file containinig gene regulatory regions", action = "store")
parser.add_argument("-t","--tf", help="please provide a bed file containing transcription factor peak regions", action = "store")
parser.add_argument("-o","--output", help="please provide an output file name", action = "store")

args = parser.parse_args()
files_path = args.output + "_results/"
tf_path = "Data/"
gene_file_nm = files_path + args.file + ".csv"
gene_df = pd.read_csv(gene_file_nm)
tf_file_nm = tf_path + args.tf

if gene_file_nm.find("promoter") != -1:
    gene_df = gene_df[["chromosome_name", "promoter_start", "promoter_end", "ensembl_gene_id"]]
    o_name = files_path + args.output + "_promoter_tfbs.csv"
else:
    gene_df = gene_df[["chromosome_name", "reg_start", "reg_end", "ensembl_gene_id"]]
    o_name  = files_path + args.output + "_tfbs.csv"

gene_df.to_csv("Gene_temp.bed", header = None, sep = "\t", index  = False)


try:
    os.system("bedtools intersect -a " + tf_file_nm + " -b Gene_temp.bed -wa -wb -f 1.0 > TF_overlap_temp.bed")
    overlap_df = pd.read_csv("TF_overlap_temp.bed", sep = "\t", header = None)
    overlap_df = overlap_df.iloc[:, [0,1,2,3,7]]
    overlap_df.columns = ["Chromosome", "TF_Start", "TF_End", "TF_Name","Ensembl"]
    overlap_df = overlap_df[["Ensembl", "Chromosome",  "TF_Start", "TF_End", "TF_Name"]]
    overlap_df.to_csv(o_name, index = False)
except:
    print("No overlapping TF peaks  were found.")

os.system("rm -rf *temp*.bed")






	



