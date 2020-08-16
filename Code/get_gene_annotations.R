library("biomaRt")
library("optparse")
library("dplyr")
option_list <- list(make_option(c("-g", "--gene"), type="character", default=NULL,
                                 help="file containing genes", metavar="character"),
                   make_option(c("-o", "--output"), type="character", default=NULL,
                                             help="Output file name", metavar="character"))

files_path <- "Data/"
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

out_path = paste(opt$output,"_results/",sep = "")

ensembl_df <- read.csv(paste(files_path,opt$gene,sep = ""), header = F)
ensembl_37 <- useEnsembl(biomart="ensembl",GRCh=37, dataset = "hsapiens_gene_ensembl")
gene_info_df <- getBM(c("ensembl_gene_id", "strand", "chromosome_name","ensembl_transcript_id","start_position","end_position","transcript_length","transcription_start_site","ensembl_exon_id", "transcript_start","transcript_end","exon_chrom_start", "exon_chrom_end"),"ensembl_gene_id", c(as.character(ensembl_df[,1])), ensembl_37)
write.csv(gene_info_df,paste(out_path, opt$output,"_gene_body.csv", sep = ""),  row.names  =F)

