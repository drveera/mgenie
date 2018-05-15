#!/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

expr.file <- args[1]
module.file <- args[2]
out.file <- args[3]

if(FALSE){
  module.file <- "chocolate.genes"
  expr.file <- "sample.raw.RDS.0.pcovar.residuals.formatted.RDS.formatted"
}

library(data.table)

expr <- fread(expr.file)
samples <- expr$sample
expr <- expr[,-1,with=FALSE]
expr <- as.matrix(expr)

modulegenes <- readLines(module.file)
expr <- expr[,intersect(colnames(expr),modulegenes)]

expr.pca <- prcomp(expr,retx = TRUE)
expr.pcs <- expr.pca$x[,1:5]
expr.pcs <- as.data.table(expr.pcs)
modulename <- gsub(".mgenes","",basename(module.file))
names(expr.pcs) <- paste0(modulename,"_",names(expr.pcs))
for(i in 1:5){
  df <- cbind(sample=samples,expr.pcs[,i,with=FALSE])
  fwrite(df,paste0(out.file,i,".expr"), na="NA")
  df.pheno <- cbind(sample=samples,df)
  fwrite(df.pheno,paste0(out.file,i,".pheno"), na="NA", sep="\t", col.names=FALSE)
}
expr.pca.importance <- data.table(t(summary(expr.pca)$importance))
fwrite(expr.pca.importance,paste0(out.file,".importance"),sep="\t",na="NA")

