#!/bin/R





library(Rtsne)
library(edgeR)
library(gplots)
library(plotrix)
library(monocle)
source("scrna_lib.r")





argv <- commandArgs(trailingOnly=T)


express_matrix <- read.csv("express_matrix.txt", header = T, sep = "\t")

groups.cl2 <- read.csv("groups.cl2", header = F, sep = "\t")[,1]

clustcol_2 <- read.csv("clustcol_2", header = F, sep = "\t")[,1]

for (i in 1: length(clustcol_2))
{
    
    group_1 <- colnames(express_matrix)[which(groups.cl2 == clustcol_2[i])] ### HFD
    
    group_2 <- colnames(express_matrix)[which(groups.cl2 != clustcol_2[i])] ### ND
    
    
    cluster_express_matrix_pre <- express_matrix[,c(group_1, group_2)]
    
    ################################ TMM normalization ###################################
    
    
    cluster_express_matrix <- tmm(cluster_express_matrix_pre)
    
    ################# monocle pipeline #################
    
    ######################## cal rpkm ########################
    
    gene_length <- read.csv("genes_length.txt", header = T, sep = "\t")
    
    d <- DGEList(counts = cluster_express_matrix)
    #cpm(d)
    #cpm(d,log=TRUE)
    
    d$genes$Length <- gene_length[,1]
    
    rpkm_express <- rpkm(d)
    
    #select_express_matrix <- cluster_express_matrix[which(rowMeans(rpkm_express) > 1), ]
    
    #select_rpkm_express_matrix <- rpkm_express[which(rowMeans(rpkm_express) > 1), ]
    
    
    rpkm_express_symbol <- rpkm_express
    
    gene_id <- read.csv("genes_title.txt", header = T, sep = "\t")
    
    rownames(rpkm_express_symbol) <- gene_id[,1]
    
    ################################# calculate Differentiall expression #######################################
    
    
    group <- c(replicate(length(group_1),"HFD"), replicate(length(group_2),"ND"))
    annot <- data.frame(group) ######## 1 vs 7 big
    
    rownames(annot) <- colnames(rpkm_express_symbol)
    
    
    colnames(annot) <- "samples"
    
    pd <- new("AnnotatedDataFrame", data = annot)
    
    
    fd_data <- data.frame(rownames(rpkm_express_symbol)) ###### remove the genes with NA symbol
    rownames(fd_data) <- rownames(rpkm_express_symbol)
    
    colnames(fd_data) <- "gene_short_name"
    
    
    fd <- new("AnnotatedDataFrame", data = fd_data)
    
    HSMM <- newCellDataSet(as.matrix(rpkm_express_symbol), phenoData = pd, featureData = fd)
    HSMM <- detectGenes(HSMM, min_expr = 0.1)
    
    # HFD_quality_num <- length(rpkm_express_symbol[which(rpkm_express_symbol[1:length(group_1)] > 0.1)])
    
    #ND_quality_num <- length(rpkm_express_symbol[which((length(group_1)+1):(length(group_1)+length(group_2)) > 0.1)])
    
    
    HSMM <- estimateSizeFactors(HSMM)
    
    diff_ND_vs_HFD <- differentialGeneTest(HSMM, fullModelFormulaStr="~samples") #########
    
    match_index_2 <- match(rownames(diff_ND_vs_HFD), rownames(rpkm_express_symbol))
    log2fc <- apply(rpkm_express_symbol, 1, function(b) {log2(mean(b[1:length(group_1)])) - log2(mean(b[(length(group_1)+1):(length(group_1)+length(group_2))]))})######## 1 vs 7 big
    
    HFD_num_qualifed <- apply(rpkm_express_symbol, 1, function(b) {length(b[which(b[1:length(group_1)] > 0.1)])})
    ND_num_qualifed <- apply(rpkm_express_symbol, 1, function(b) {length(b[which(b[(length(group_1)+1):(length(group_1)+length(group_2))] > 0.1)])})
    
    result_matrix <- data.frame(diff_ND_vs_HFD,HFD_num_qualifed[match_index_2], ND_num_qualifed[match_index_2], as.matrix(log2fc)[match_index_2,], rpkm_express_symbol[match_index_2,]) #############
    
    ################################# converting gene id #######################################
    
    # annotate the genes
    keys <- as.character(rownames(result_matrix))
    
    if (argv[1] == "Human")
	{
		library(org.Hs.eg.db)
		annotations <- select(org.Hs.eg.db, keys,  c("SYMBOL","GENENAME"), keytype="ENSEMBL")
	}

	if (argv[1] == "Mouse")
	{
		library(org.Mm.eg.db)
		annotations <- select(org.Mm.eg.db, keys,  c("SYMBOL","GENENAME"), keytype="ENSEMBL")
	}

	if (argv[1] == "Rat")
	{
		library(org.Rn.eg.db)
		annotations <- select(org.Rn.eg.db, keys,  c("SYMBOL","GENENAME"), keytype="ENSEMBL")
	}


  
    match_index <- match(annotations[,1], keys)
    output <- data.frame(annotations, result_matrix[match_index,])
    
    filename = paste(clustcol_2[i], "_vs_rest.txt", sep="")
    write.table(output, filename, row.names = T, col.names=T, quote = F,sep="\t")  #############
    system(paste("sed 's/\\./,/g' ", clustcol_2[i], "_vs_rest.txt > ", clustcol_2[i], "_vs_rest_comma.txt", sep = ""))
    
}
