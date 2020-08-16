library(pandaR)
library(optparse)

option_list = list( make_option(c("-m", "--motif"), type="character", default=NULL,
			                     help="motif file name in csv format", metavar="character"),
		  	make_option(c("-p", "--ppi"), type="character", default=NULL,
				                    help="ppi file name in csv format", metavar="character"),
		  	make_option(c("-e", "--expression"), type="character", default=NULL,
				                    help="expression file name in csv format", metavar="character"),
		  	make_option(c("-d", "--hammingdistance"), type="numeric", default=0.001,
				                    help="hamming distance [default= %default]", metavar="character"),
		  	make_option(c("-o", "--output"), type="character", default=NULL,
				                    help="output file prefix", metavar="character"))
				            


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
files_path <- "Data/"
out_path <-  paste(opt$output,"_results/",sep = "")
mtf <- read.csv(paste(out_path,opt$motif,sep = ""))
colnames(mtf) <- c("tf_name", "hgnc_symbol", "Adj")

if (opt$ppi == "none"){
	tf_names <- unique(mtf$tf_name)
	ppi <- diag(length(tf_names))
	rownames(ppi) <- tf_names;colnames(ppi) <- tf_names
}else{
	ppi <- read.csv(paste(files_path,opt$ppi,sep = ""))
}

if (opt$expression == "none"){
	gene_list = unique(mtf$hgnc_symbol)
	expr2 <- as.data.frame(cbind(c(rep(1,length(gene_list))), c(rep(1, length(gene_list)))))
	row.names(expr2) <- gene_list
}else{
	expr <-  read.csv(paste(files_path,opt$expression,sep = ""))
	expr2 <- expr[,2:ncol(expr)]
	row.names(expr2) <- c(as.character(expr[,1]))
}

hp <- panda(mtf, expr2,ppi,hamming = opt$hammingdistance,progress = T,remove.missing.ppi = T, remove.missing.motif = T, remove.missing.genes = T)
o_name = paste(out_path,opt$output,"_reg_net.tsv",sep="")
write.table(hp@regNet, file = o_name, sep = "\t")



