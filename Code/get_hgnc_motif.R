library("biomaRt")
library("optparse")
library("dplyr")
option_list <- list(make_option(c("-t", "--tfbs"), type="character", default=NULL,
				 help="TFBS file name in csv format", metavar="character"),
		   make_option(c("-o", "--output"), type="character", default=NULL,
			                     help="output file name prefix", metavar="character"),
                   make_option(c("-y" , "--type"), action = "store", help = "provide option for the type of input tfbs (weighted or unweighted)",default = NULL, metavar = "character"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

files_path <- paste(opt$output,"_results/",sep = "")
out_path = files_path

motif_df <- read.csv(paste(files_path,opt$tfbs,sep = ""))

ensembl_37 <- useEnsembl(biomart="ensembl",GRCh=37, dataset = "hsapiens_gene_ensembl")

if (opt$type == "weighted") {
	ensembl_df <- motif_df
	colnames(ensembl_df) <- c("TF_Name", "ensembl_gene_id","Adj")
	hg_df <- getBM(c("ensembl_gene_id", "hgnc_symbol"),"ensembl_gene_id", unique(ensembl_df$ensembl_gene_id), ensembl_37)
	ensembl_merge_df <- inner_join(ensembl_df, hg_df)
	
}else{
	ensembl_df <- distinct(motif_df[,c("TF_Name", "Ensembl")])
	colnames(ensembl_df) <- c("TF_Name", "ensembl_gene_id")
	#distinct(motif_df[,c("TF_Name", "Ensembl")])
	hg_df <- getBM(c("ensembl_gene_id", "hgnc_symbol"),"ensembl_gene_id", unique(ensembl_df$ensembl_gene_id), ensembl_37)
	ensembl_merge_df <- inner_join(ensembl_df, hg_df)
	ensembl_merge_df$Adj = c(rep(1,nrow(ensembl_merge_df)))
}

ensembl_merge_df <- ensembl_merge_df[,c("TF_Name", "hgnc_symbol","Adj")]

write.csv(ensembl_merge_df,paste(files_path,opt$output,"_motif.csv", sep = ""), row.names  =F)

