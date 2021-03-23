import scipy
import pandas as pd
import numpy as np
import glob
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
import itertools as itr
import multiprocessing as mp
import random
from sklearn.linear_model import ElasticNetCV,RidgeCV
from sklearn.metrics import mean_squared_error
import multiprocessing as mp
from functools import partial
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", help="please provide a file name containing input file names for running the models", action = "store")
parser.add_argument("-o","--output", help="please provide an output file name", action = "store")
parser.add_argument("-e","--expression", help="please provide the file name containing expression fpkm values", action = "store")
parser.add_argument("-r","--random", help="please provide the number of random instances for which the ENET models will be run",type  = int, nargs = "?", const = 20,  action = "store")
parser.add_argument("-m","--multiprocess", help = "select this flag along with all flag to do multiprocessing", action = "store_true")
parser.add_argument("-n", "--cpu",type = int, help = "number of cpu cores to use", action = "store")
args = parser.parse_args()

input_file_path  = "Data/"
expr_df  =pd.read_csv(input_file_path + args.expression,header = None)
out_path =  "pred_results/"

if not os.path.exists(out_path):
	os.mkdir(out_path)


def get_genes(y):
        panda_scores_df =pd.read_csv(y, sep = "\t")
        panda_scores_t = panda_scores_df.T
        panda_scores_t.insert(0,"hgnc_symbol", list(panda_scores_t.index))
        panda_scores_t1 = panda_scores_t
        return(set(panda_scores_t1.hgnc_symbol))

file_name = args.input

with open(file_name) as fl:
    file_list  = fl.read().splitlines()
fl.close()

file_list = [i for i in file_list if len(i) > 1]

print("Finding common genes")

if args.multiprocess:
	pool = mp.Pool(args.cpu)
	list_genes = list(pool.map(get_genes, file_list))
	del(pool)
else:
        list_genes = list(map(get_genes, file_list))

common_genes = set(list_genes[0])
for i in list_genes[1:]:
        common_genes.intersection_update(i)

expr_df1 = expr_df.loc[expr_df[0].isin(list(common_genes))]
expr_df1.columns = ["hgnc_symbol", "fpkm"]
common_genes_f = list(sorted(expr_df1.hgnc_symbol))
print(len(common_genes_f), " common genes were found.")

def get_preds(y,x):
    k1 = y.split("/")
    len_k = len(k1) - 1
    if y.find("reg") != -1:
        panda_scores_df = pd.read_csv(y, sep = "\t")
        panda_scores_t = panda_scores_df.T
        panda_scores_t.insert(0,"hgnc_symbol", list(panda_scores_t.index))
        #panda_scores_t = panda_scores_t.loc[panda_scores_t.hgnc_symbol.isin(common_genes_f)]
        k2 = "_".join(k1[len_k].split("_")[:3])
    else:
        panda_scores_df = pd.read_csv(y)
        panda_scores_t = panda_scores_df    
        k2 = "_".join(k1[len_k].split("_")[:2])
    panda_scores_t = panda_scores_t.loc[panda_scores_t.hgnc_symbol.isin(common_genes_f)]
    mrg_df = expr_df1.merge(panda_scores_t,on = "hgnc_symbol", how = "inner")
    np_scores = np.array(mrg_df.iloc[:,2:])
    np_scores_scaled  = MinMaxScaler().fit_transform(np_scores)
    np_expr = np.array(mrg_df.iloc[:,1])
    np_expr_scaled = StandardScaler().fit_transform(np_expr.reshape((-1,1)))
    np_expr_scaled = np_expr_scaled.reshape((-1,))
    x_train,x_test,y_train,y_test = train_test_split(np_scores_scaled,np_expr_scaled, test_size = 0.20,random_state = x)
    enet_model = ElasticNetCV(cv =20, max_iter = 10000, tol = 0.001)
    enet_model.fit(x_train,y_train)
    pred_enet = enet_model.predict(x_test)       
    corr_enet = np.corrcoef(y_test,pred_enet.reshape(-1,))[1,0]
    mse_enet = mean_squared_error(y_test,pred_enet)
    mrg_list = [x,k2,corr_enet,mse_enet]
    return(mrg_list)


random_states_f = random.sample(range(1000,10000000), args.random)

o_name = out_path + args.output + "_enet_pred_eval.csv"
pool = mp.Pool(args.cpu)
file_name = args.input

with open(file_name) as fl:
    file_list  = fl.read().splitlines()
fl.close()

fl_name = file_list[0]
if fl_name.find("reg") != -1:
        panda_df = pd.read_csv(fl_name, sep = "\t")
        panda_t = panda_df.T
        panda_t.insert(0,"hgnc_symbol", list(panda_t.index))
else:
        panda_df = pd.read_csv(fl_name)
        panda_t = panda_df   
panda_t = panda_t.loc[panda_t.hgnc_symbol.isin(common_genes_f)]
gene_train, gene_test = train_test_split(list(panda_t.hgnc_symbol), test_size = 0.20, random_state = 100)
print("Running predictions for ", str(args.random),"iterations.")
print("Building the models using ", str(len(gene_train)), "genes and testing it on ", str(len(gene_test)))

def get_pred_state(n):
    partial_run_pred  =  partial(get_preds,x = n)
    df_list = list(pool.map(partial_run_pred, file_list))
    df_f = pd.DataFrame(df_list)
    df_f.columns =  ["Random_state", "Input_matrix", "Corr_ENET", "MSE_ENET"]
    return(df_f)

if args.multiprocess:
    dlist = list(map(get_pred_state, random_states_f))
    dl = pd.concat(dlist)
    
else:
    dfl = []
    for i in random_states_f:
        for k in file_list:
            dfl.append(get_preds(i,k))
    dl = pd.DataFrame(df_list)
    dl.columns =  ["Random_state", "Input_matrix", "Corr_ENET", "MSE_ENET"]
    
dl.to_csv(o_name, index= False)














