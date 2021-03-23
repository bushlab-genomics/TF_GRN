import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-p","--promoter", help="please provide a file containing promter tfbs", action = "store")
parser.add_argument("-t","--tf", help="please provide filename containing all tfbs", action = "store")
parser.add_argument("-o","--output", help="please provide an output file name", action = "store")

args = parser.parse_args()
out_path = args.output + "_results/"
all_tfbs = pd.read_csv(out_path + args.tf)
prom_tfbs = pd.read_csv(out_path + args.promoter)
prom_tfbs = prom_tfbs.drop_duplicates()
all_tfbs = all_tfbs.drop_duplicates()
common_genes = list(set.intersection(set(all_tfbs.Ensembl), set(prom_tfbs.Ensembl)))
all_tfbs = all_tfbs.loc[all_tfbs.Ensembl.isin(common_genes)]
prom_tfbs = prom_tfbs.loc[prom_tfbs.Ensembl.isin(common_genes)]

all_tfbs["reg_index"] = list(all_tfbs.index)
prom_tfbs_ind = prom_tfbs.merge(all_tfbs, on = list(prom_tfbs), how = "inner")
distal_tfbs = all_tfbs.loc[~all_tfbs.reg_index.isin(list(prom_tfbs_ind.reg_index))]
distal_tfbs = distal_tfbs.iloc[:,:5]
o_name = out_path + args.output + "_distal_tfbs.csv"

distal_tfbs.to_csv(o_name, index= False)


