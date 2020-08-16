import pandas as pd
import numpy as np
import argparse
import multiprocessing as mp
import os

parser = argparse.ArgumentParser()
parser.add_argument("-i","--inputs", help="please provide an input filenames for input gene list, peak file,ctcf file, panda ppi and panda expression",nargs = 5, action = "store")
parser.add_argument("-t", "--type", help = "please provide the type of grn to be contructed [intron, hic, weighted, promoter, distal, grn]",action =  "store")
parser.add_argument("-w", "--window", help = "please provide the cis regulatory window(bp) for finding regulatory regions(default = 50000)", nargs = "?", const = 50000,type = int, action = "store")
parser.add_argument("-p", "--promoterwindow", help = "please provide the cis regulatory window(bp) for finding  promoter regions(default = 5000)",nargs = "?", const = 5000,  type = int, action = "store")
parser.add_argument("-s","--scalingfactor",  help="please provide a scaling factor in the range [0,1] for HiC normalized motif(default = 1.0)",type = float, const = 1.0, nargs = "?", action = "store")
parser.add_argument("-u","--upweighting",  help="select this flag if you want  to upweight the promoter TFBS in the HiC analysis", action = "store_true")
parser.add_argument("-c","--contacts", help = "please provide the name of the file containing hic contacts", action = "store")
parser.add_argument("-d", "--hammingdistance", help = "please provide the hamming distance for PANDA GRN convergence(default = 0.001)", nargs = "?", const = 0.001, type = float, action = "store")
parser.add_argument("-m","--multiprocessing", help = "select this flag along with all flag to do multiprocessing", action = "store_true")
parser.add_argument("-n", "--cpu",type = int, help = "number of cpu cores to use", action = "store")
parser.add_argument("-o","--output", help="please provide an output file prefix", action = "store")
args = parser.parse_args()
path =  args.output + "_results"
if not os.path.exists(path):
	os.mkdir(path)

gene_list = args.inputs[0]
peak_file = args.inputs[1]
ctcf_file = args.inputs[2]
panda_ppi  = args.inputs[3]
panda_expression = args.inputs[4]

print("Generating gene annotations.")

os.system("Rscript Code/get_gene_annotations.R -g " + gene_list + " -o " + args.output)

if args.multiprocessing:
	if args.type == "promoter":
		print("Defining promoter regions.")
		os.system("python Code/get_promoter_regions.py -i " + args.output + "_gene_body -m -n " + str(args.cpu) + " -w " + str(args.promoterwindow) + " -o " + args.output)
		print("Finding promoter TFBS.")
		os.system("python Code/get_tfbs.py -f " + args.output + "_promoter_regions -t " + peak_file + " -o " + args.output)
		print("Generating promoter motif network.")
		os.system("Rscript Code/get_hgnc_motif.R -t " + args.output + "_promoter_tfbs.csv -y unweighted -o " + args.output)
	else:
		print("Defining cis regulatory regions.")
		os.system("python Code/get_cis_regulatory_regions.py -i " + args.output + "_gene_body -m -n " +  str(args.cpu) + " -c " + ctcf_file + " -w " + str(args.window) + " -o " + args.output)
		print("Finding TFBS.")
		os.system("python Code/get_tfbs.py -f " + args.output + "_reg_regions -t " + peak_file +  " -o " + args.output)
		if args.type == "distal":
			print("Defining promoter regions.")
			os.system("python Code/get_promoter_regions.py -i " + args.output + "_gene_body -m -n " + str(args.cpu) + " -w " + str(args.promoterwindow) + " -o " + args.output)	
			print("Finding promoter TFBS.")
			os.system("python Code/get_tfbs.py -f " + args.output +"_promoter_regions -t " + peak_file + " -o " + args.output) 
			print("Finding distal TFBS.")
			os.system("python Code/get_distal_tfbs.py -p " + args.output + "_promoter_tfbs.csv -t " + args.output + "_tfbs.csv -o " + args.output)
			print("Generating distal motif network.")
			os.system("Rscript Code/get_hgnc_motif.R -t " + args.output + "_distal_tfbs.csv -y unweighted -o " + args.output)
		elif args.type == "intron":
			print("Defining intron regions.")
			os.system("python Code/get_intron_regions.py -i " + args.output + "_gene_body.csv -r " + args.output + "_reg_regions.csv -m -n " + str(args.cpu) + " -o " + args.output)
			print("Finding intronic TFBS.")
			os.system("python Code/get_tfbs.py -f " + args.output + "_intron_regulatory_regions -t " + peak_file + " -o " + args.output)
			print("Generating intron motif network.")
			os.system("Rscript Code/get_hgnc_motif.R -t " + args.output + "_tfbs.csv -y unweighted -o " + args.output)
		elif args.type == "hic":
			print("Defining promoter regions.")
			os.system("python Code/get_promoter_regions.py -i " + args.output + "_gene_body -m -n " + str(args.cpu) + " -w " + str(args.promoterwindow) + " -o " + args.output)
			print("Finding promoter TFBS.")
			os.system("python Code/get_tfbs.py -f " + args.output + "_promoter_regions -t " + peak_file +  " -o " + args.output)			
			print("Overlapping Hi-C contact matrix with the TFBS and promoter regions")
			os.system("python Code/hic_join_bedtools.py -p " + args.output + "_promoter_regions.csv -c " + args.contacts + " -t " + args.output + "_tfbs.csv -o " + args.output)
			print("Normalizing Hi-C motif network")
			if args.upweighting:
				os.system("python Code/get_scaled_contacts.py -t " + args.output + "_tfbs.csv -c " + args.output + "_hic_merge.csv -u -p " + args.output + "_promoter_tfbs.csv -s " + str(args.scalingfactor) + " -o " + args.output)   
			else:
				os.system("python Code/get_scaled_contacts.py -t " + args.output + "_tfbs.csv -c " + args.output + "_hic_merge.csv -p " + args.output + "_promoter_tfbs.csv -s " + str(args.scalingfactor) + " -o " + args.output)
			os.system("Rscript Code/get_hgnc_motif.R -t " + args.output + "_hic_tfbs.csv -y weighted -o " + args.output)
		elif args.type == "weighted":
			print("Generating TF-TG weights.")
			os.system("python Code/get_tfbs_weights.py -i " + args.output + "_tfbs.csv -m -n " + str(args.cpu) + " -o " + args.output)
			print("Generating weighted motif network")
			os.system("Rscript Code/get_hgnc_motif.R -t " + args.output + "_weighted_tfbs.csv -y weighted -o " + args.output)
		elif args.type == "grn":
			print("Generating motif network.")
			os.system("Rscript Code/get_hgnc_motif.R -t " + args.output + "_tfbs.csv -y unweighted -o " + args.output)

else:
	if args.type == "promoter":
		print("Defining promoter regions.")
		os.system("python Code/get_promoter_regions.py -i " + args.output + "_gene_body -w " + str(args.promoterwindow) + " -o " + args.output)
		print("Finding promoter TFBS.")
		os.system("python Code/get_tfbs.py -f " + args.output + "_promoter_regions -o " + args.output)
		print("Generating promoter motif network.")
		os.system("Rscript Code/get_hgnc_motif.R -t " + args.output + "_promoter_tfbs.csv -y unweighted -o " + args.output)
	else:
		print("Defining cis regulatory regions.")
		os.system("python Code/get_cis_regulatory_regions.py -i " + args.output + "_gene_body -c " + ctcf_file + " -w " + str(args.window) + " -o " + args.output)
		print("Finding TFBS.")
		os.system("python Code/get_tfbs.py -f " + args.output + "_reg_regions -t " + peak_file +  " -o " + args.output)
		if args.type == "distal":
			print("Defining promoter regions.")
			os.system("python Code/get_promoter_regions.py -i " + args.output + "_gene_body -w " + str(args.promoterwindow) + " -o " + args.output)
			print("Finding promoter TFBS.")
			os.system("python Code/get_tfbs.py -f " + args.output + "_promoter_regions -o " + args.output)
			print("Finding distal TFBS.")
			os.system("python Code/get_distal_tfbs.py -p " + args.output + "_promoter_tfbs.csv -t " + args.output + "_tfbs.csv -o " + args.output)
			print("Generating distal motif network.")
			os.system("Rscript Code/get_hgnc_motif.R -t " + args.output + "_distal_tfbs.csv -y unweighted -o " + args.output)
                        
		elif args.type == "intron":
			print("Defining intron regions.")
			os.system("python Code/get_intron_regions.py -i " + args.output + "_gene_body.csv -r " + args.output + "_reg_regions.csv -o " + args.output)
			print("Finding intronic TFBS.")
			os.system("python Code/get_tfbs.py -f " + args.output + "_intron_regulatory_regions -t " + peak_file + " -o " + args.output)
			print("Generating intron motif network.")
			os.system("Rscript Code/get_hgnc_motif.R -t " + args.output + "_tfbs.csv -y unweighted -o " + args.output)
		elif args.type == "hic":
			print("Defining promoter regions.")
			os.system("python Code/get_promoter_regions.py -i " + args.output + "_gene_body -w " + str(args.promoterwindow) + " -o " + args.output)
			print("Finding promoter TFBS.")
			os.system("python Code/get_tfbs.py -f " + args.output + "_promoter_regions -t " + peak_file +  " -o " + args.output)
			print("Overlapping Hi-C contact matrix with the TFBS and promoter regions")
			os.system("python Code/hic_join_bedtools.py -p " + args.output + "_promoter_regions.csv -c " + args.contacts + " -t " + args.output + "_tfbs.csv -o " + args.output)
			print("Normalizing Hi-C motif network")
			if args.upweighting:
				os.system("python Code/get_scaled_contacts.py -t " + args.output + "_tfbs.csv -c " + args.output + "_hic_merge.csv -u -p " + args.output + "_promoter_tfbs.csv -s " + str(args.scalingfactor) + " -o " + args.output)
			else:
				os.system("python Code/get_scaled_contacts.py -t " + args.output + "_tfbs.csv -c " + args.output + "_hic_merge.csv -p " + args.output + "_promoter_tfbs.csv -s " + str(args.scalingfactor) + " -o " + args.output)
			os.system("Rscript Code/get_hgnc_motif.R -t " + args.output + "_hic_tfbs.csv -y weighted -o " + args.output)
		elif args.type == "weighted":
			print("Generating TF-TG weights.")
			os.system("python Code/get_tfbs_weights.py -i " + args.output + "_tfbs.csv -o " + args.output)
			print("Generating weighted motif network")
			os.system("Rscript Code/get_hgnc_motif.R -t " + args.output + "_weighted_tfbs.csv -y weighted -o " + args.output)
		elif args.type == "grn":
			print("Generating motif network.")
			os.system("Rscript Code/get_hgnc_motif.R -t " + args.output + "_tfbs.csv -y unweighted -o " + args.output)

print("Generating PANDA GRN.")
os.system("Rscript Code/run_panda.R -m " + args.output + "_motif.csv" + " -p " + panda_ppi + " -e " + panda_expression + " -d " + str(args.hammingdistance) +  " -o " + args.output)

