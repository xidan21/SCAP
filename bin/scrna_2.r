#!/bin/R




library(Rtsne)
library(edgeR)
library(gplots)
library(plotrix)
library(monocle)
source("scrna_lib.r")




argv <- commandArgs(trailingOnly=T)

NORM <- read.csv("../source/NORM.txt", header = T, sep = "\t")
RPKM <- read.csv("../source/RPKM.txt", header = T, sep = "\t")

ercc <- grep("ERCC",rownames(RPKM))
nD.ercc <- apply(RPKM[ercc,],2,detect.genes)
nD.tx <- apply(RPKM[-ercc,],2,detect.genes) ## rpkm > 1
ercc.ratio<-colSums(RPKM[-ercc,])/colSums(RPKM[ercc,])


passed<-which(nD.tx >=as.numeric(argv[2]) & nD.tx<=as.numeric(argv[3]))


length(passed)

nS <- ncol(RPKM) # number of samples

length(passed)/nS

NORM <- NORM[,passed]
RPKM <- RPKM[,passed]
nG<-nrow(NORM)
nS<-ncol(NORM)

##################### top variable genes ####################

if (argv[6] == "all")
{
	
	NORM <- NORM
} else {

	NORM <- top_var(NORM, as.numeric(argv[6]))
}

##################### plot PCA ####################



nD.tx <- apply(NORM,2,detect.genes) # with RPKM > 1 and remove ercc and remove the bad cells
cc <- colorRampPalette(c("yellow","red","black"))
col.nDet <- convert.to.color(nD.tx,cc)


PC <- run.pca(NORM)
pca.plot(PC,pch=16,main="PCA, color by detected genes",col=col.nDet$cols)


################## plot tSNE #####################

tsne2 <- Rtsne(t(log2(NORM+1)), perplexity = as.numeric(argv[5]))

par(mfrow=c(2,2),mar=c(2,2,2,2),oma=c(1,1,1,1))

plot(tsne2$Y,col=col.nDet$cols,pch=16)

# test hclust on the tSNE matrix
hc2 <- hclust(dist(tsne2$Y),method="ward.D2")
groups.cl <- cutree(hc2,as.numeric(argv[4]))

hc2_cluster <- rect.hclust(hc2, k=as.numeric(argv[4])) # define groups

clustcol_2 <- palette()[1:length(hc2_cluster)]

groups.cl2 <- mgsub(as.character(seq(as.numeric(argv[4]))), clustcol_2, groups.cl) 

colplot.clust(hc2,lab.col=groups.cl2)

#dev.off()

pdf (file = paste("../result/tsne_plot_", argv[5],".pdf", sep=""))

plot(tsne2$Y,col=groups.cl2,pch=16, main = paste("perplexity", argv[5]))

dev.off()
#text(tsne2$Y, labels = colnames(RPKM), cex = 0.3)

write.table(tsne2$Y, "../result/tsne2_y.txt", col.names = T, row.names = T, quote = F, sep = "\t")


write.table(groups.cl2,"../result/groups.cl2.txt", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(RPKM,"../result/RPKM.txt", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(clustcol_2,"../source/clustcol_2", col.names = F, row.names = F, quote = F, sep = "\t")

################## plot single gene tSNE #####################


keys <- as.character(rownames(RPKM))

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
RPKM <- RPKM[match_index,]

rownames(RPKM) <- make.names(annotations[,2], unique = T)

convert.to.color <- function(x,colscale,col.range = NULL){
  x.range <- range(na.omit(x))
  by=0.1
  if (is.null(col.range)){
    col.range <- seq(x.range[1],x.range[2],by=by)
  }
  col.def <- colscale(length(col.range))
  col.idx <- round((x-x.range[1])/by)+1
  col.idx[col.idx > length(col.range)] <- length(col.range)
  cols <- col.def[col.idx]
  return(list(cols=cols,col.def=col.def,col.range=col.range))
}

candidate_genes <- read.csv("../source/marker_genes.txt", header = F)[,1]

print(candidate_genes)

for (i in 1:length(candidate_genes))
{

	filename = paste("../result/",candidate_genes[i], '_',argv[5],'_marker_genes.pdf', sep = "") ###################### need change #########################

	pdf (file= filename)
	bla.cl <- convert.to.color(unlist(RPKM[as.character(candidate_genes[i]),]),cc)

	print (i)
	print (candidate_genes[i])
	plot(tsne2$Y,col=bla.cl$cols,pch=16, main = paste(candidate_genes[i], "perplexity", argv[5])) ###################### need change #########################
	
	dev.off()

}

