#!/bin/R



library(Rtsne)
library(DESeq2)
library(edgeR)
library(gplots)
library(plotrix)
#library(monocle)
source("scrna_lib.r")




argv <- commandArgs(trailingOnly=T)

print(argv[1])
# argv[1] = "express_matrix.txt"
# argv[2] = "genes_title.txt"
# argv[3] = "genes_length.txt"


express_matrix <- read.csv(argv[1], header = T, sep = "\t")
gene_id <- read.csv(argv[2], header = T, sep = "\t")
rownames(express_matrix) <- gene_id[,1]

######################## deseq2 normalization ########################

normalized_matrix <- deseq2_normalizaton(express_matrix, 2, ncol(express_matrix)-2)


######################## cal rpkm ########################

gene_length <- read.csv(argv[3], header = T, sep = "\t")

d <- DGEList(counts=express_matrix)
#cpm(d)
#cpm(d,log=TRUE)

d$genes$Length <- gene_length[,1]

rpkm_express_matrix <- rpkm(d)
gene_id[,1] -> gene.name


RPKM <- rpkm_express_matrix

NORM <- normalized_matrix

rownames(RPKM)<- gene.name

rownames(NORM)<- gene.name


######################## plot ercc ########################

ercc <- grep("ERCC",rownames(RPKM))
nD.ercc <- apply(RPKM[ercc,],2,detect.genes)
nD.tx <- apply(RPKM[-ercc,],2,detect.genes) ## rpkm > 1
ercc.ratio<-colSums(RPKM[-ercc,])/colSums(RPKM[ercc,])

pdf (file="../result/detected_transcripts.pdf")
par(mfrow=c(2,1),mar=c(3,3,1,1),cex=0.6)
hist(nD.tx,n=100,main="Detected genes")
hist(nD.ercc,n=100,main="Detected ERCC")



dev.off()

write.table(RPKM,"../source/RPKM.txt", col.names = T, row.names = T, quote = F, sep = "\t")
write.table(NORM,"../source/NORM.txt", col.names = T, row.names = T, quote = F, sep = "\t")
