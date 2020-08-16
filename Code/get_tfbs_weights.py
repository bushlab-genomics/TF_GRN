import pandas as pd
import numpy as np
import argparse
import multiprocessing as mp

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", help="please provide a file name containing input file name containing the TFBS", action = "store")
parser.add_argument("-o","--output", help="please provide the output file name", action = "store")
parser.add_argument("-m","--multiprocessing", help = "select this flag along with all flag to do multiprocessing", action = "store_true")
parser.add_argument("-n", "--cpu",type = int, help = "number of cpu cores to use", action = "store")
args = parser.parse_args()
files_path = args.output + "_results/"

tfbs_df = pd.read_csv(files_path + args.input)

def get_normalized_motif(x):
    gn_df = tfbs_df.loc[tfbs_df.Ensembl == x]
    gn_df_length = len(gn_df)
    gn_tf_list = list(sorted(set(gn_df.TF_Name)))
    motif_score_list = []
    for i in gn_tf_list:
        motif_score = len(gn_df.loc[gn_df.TF_Name == i])/gn_df_length
        motif_score_list.append(motif_score)
    gn_motif_df = pd.DataFrame({"TF_Name":gn_tf_list,"Ensembl":[x]*len(gn_tf_list),"Interaction_Weight":motif_score_list})
    return(gn_motif_df)


if args.multiprocessing:
	pool = mp.Pool(args.cpu)
	gene_list = sorted(set(tfbs_df.Ensembl))
	df_list = pd.concat(list(pool.map(get_normalized_motif, gene_list)))
else:
	gene_list = sorted(set(tfbs_df.Ensembl))
	df_list = pd.concat(list(map(get_normalized_motif, gene_list)))

o_name = files_path + args.output + "_weighted_tfbs.csv"
df_list.to_csv(o_name, index = False)




