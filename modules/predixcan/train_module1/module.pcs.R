#!/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

module.file <- args[1]
out.file <- args[2]
##expr.file <- gsub("_color_.*$","",module.file)
expr.file <- args[3]

if(FALSE){
  module.file <- "modules.names/Gx.WB_yellow"
  out.file <- "test.delete"
  expr.file <- "gtex.whole.blood.IDd.RDS.0.pcovar.residuals"
}

library(data.table)

expr <- t(readRDS(expr.file))
expr <- cbind(sample=rownames(expr),data.table(expr))
##expr <- fread(expr.file)
samples <- expr$sample
expr <- expr[,-c("sample"),with=FALSE]
expr <- as.matrix(expr)
rownames(expr) <- samples
##sample <- rownames(expr)
##sample <- paste0(sample,"_",sample) ####temp
##rownames(expr) <- sample
##expr <- data.table(expr.rds)
##expr <- cbind(sample=sample,expr)
##samples <- expr$sample
##expr <- expr[,-1,with=FALSE]
##expr <- as.matrix(expr)

modulegenes <- readLines(module.file)
modulegenes <- gsub("\\..*$","",modulegenes)
expr <- expr[,intersect(colnames(expr),modulegenes)]

expr.pca <- prcomp(expr,retx = TRUE)
##all PCs>1% variance
expr.pca.importance.0 <- data.table(t(summary(expr.pca)$importance))
names(expr.pca.importance.0) <- c("SD","V","PropV")
expr.pca.importance <- expr.pca.importance.0[V>=0.01]
npcs <- nrow(expr.pca.importance)
if(npcs==0){
  expr.pca.importance <- expr.pca.importance.0[1,]
}
fwrite(expr.pca.importance,paste0(out.file,".importance"),sep="\t",na="NA")
###
expr.pcs <- expr.pca$x[,1:npcs]
pc.sample <- rownames(expr.pcs)
expr.pcs <- as.data.table(expr.pcs)
##modulename <- gsub(".mgenes","",basename(module.file))
modulename <- basename(module.file)
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
if(npcs==0){
  df <- cbind(sample=pc.sample,expr.pcs[,1,with=FALSE])
  fwrite(df,paste0(out.file,0,".expr"), na="NA")
  df.pheno <- cbind(sample=pc.sample,df)
  fwrite(df.pheno,paste0(out.file,0,".pheno"), na="NA", sep="\t", col.names=FALSE)
  siggenes <- findgenes(expr,df)
  writeLines(siggenes,paste0(out.file,0,".sigGenes"))
} else {
  for(i in 1:npcs){
    df <- cbind(sample=pc.sample,expr.pcs[,i,with=FALSE])
    fwrite(df,paste0(out.file,i,".expr"), na="NA")
    df.pheno <- cbind(sample=pc.sample,df)
    fwrite(df.pheno,paste0(out.file,i,".pheno"), na="NA", sep="\t", col.names=FALSE)
    siggenes <- findgenes(expr,df)
    writeLines(siggenes,paste0(out.file,i,".sigGenes"))
  }
}




