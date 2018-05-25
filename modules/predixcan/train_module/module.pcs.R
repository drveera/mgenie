#!/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

module.file <- args[1]
out.file <- args[2]
expr.file <- gsub("_color_.*$","",module.file)

if(FALSE){
  module.file <- "sample.residuals_color_.genes"
}

library(data.table)

expr <- t(readRDS(expr.file))
sample <- rownames(expr)
sample <- paste0(sample,"_",sample) ####temp
rownames(expr) <- sample
##expr <- data.table(expr.rds)
##expr <- cbind(sample=sample,expr)
##samples <- expr$sample
##expr <- expr[,-1,with=FALSE]
##expr <- as.matrix(expr)

modulegenes <- readLines(module.file)
expr <- expr[,intersect(colnames(expr),modulegenes)]

expr.pca <- prcomp(expr,retx = TRUE)
expr.pcs <- expr.pca$x[,1:5]
pc.sample <- rownames(expr.pcs)
expr.pcs <- as.data.table(expr.pcs)
modulename <- gsub(".mgenes","",basename(module.file))
names(expr.pcs) <- paste0(modulename,"_",names(expr.pcs))

###function to find correlated genes
findgenes <- function(expr,pc.dfm){
  sample <- rownames(expr)
  expr <- data.table(expr)
  genes <- colnames(expr)
  expr <- cbind(sample=sample,expr)
  names(pc.dfm)[2] <- "pc"
  dfm <- merge(expr,pc.dfm, by="sample")
  pc <- dfm$pc
  dfm.sub <- dfm[,genes,with=FALSE]
  gcor <- apply(dfm.sub,2, function(x) return(cor.test(x,pc)$p.value))
  ##gcor <- gcor * length(gcor)
  gcor <- p.adjust(gcor,method="BH")
  siggenes <- names(gcor[gcor<0.05])
  return(siggenes)
}

for(i in 1:5){
  df <- cbind(sample=pc.sample,expr.pcs[,i,with=FALSE])
  fwrite(df,paste0(out.file,i,".expr"), na="NA")
  df.pheno <- cbind(sample=pc.sample,df)
  fwrite(df.pheno,paste0(out.file,i,".pheno"), na="NA", sep="\t", col.names=FALSE)
  siggenes <- findgenes(expr,df)
  writeLines(siggenes,paste0(out.file,i,".sigGenes"))
}
expr.pca.importance <- data.table(t(summary(expr.pca)$importance))
fwrite(expr.pca.importance,paste0(out.file,".importance"),sep="\t",na="NA")

