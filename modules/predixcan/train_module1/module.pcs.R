#!/bin/Rscript

### * arguments
args <- commandArgs(trailingOnly = TRUE)
module.file <- args[1]
out.file <- args[2]
##expr.file <- gsub("_color_.*$","",module.file)
expr.file <- args[3]

### * debug arguments
if(FALSE){
  module.file <- "modules.names/Gx.WB_yellow"
  out.file <- "test.delete"
  expr.file <- "gtex.whole.blood.IDd.RDS.0.pcovar.residuals"
}

### * libraries
library(data.table)

### * function
### ** find pc
findpc <- function(n,pcnum){
  ##shuffle rows
  expr.x <- expr.full[sample(nrow(expr.full)),]
  rownames(expr.x) <- rownames(expr.full)
  ##shuffle columns
  expr.y <- expr.full[,sample(ncol(expr.full))]
  rownames(expr.y) <- rownames(expr.full)
  ##shuffle both rows and columns
  library(picante)
  expr.xy <- randomizeMatrix(expr.full,null.model="richness",iterations=1000)
  rownames(expr.xy) <- rownames(expr.full)
  ##make a list
  expr.list <- list(expr.full,expr.x,expr.y,expr.xy)
  names(expr.list) <- c("NULL","NULLX","NULLY","NULLXY")
  ## get pcs from all
  x.pcs.l <- list()
  for(i in names(expr.list)){
    ##choose genes randomly from the columns
    x <- expr.list[[i]][,sample(1:ncol(expr.list[[i]]),n)]
    x.pca <- prcomp(x, retx = TRUE)
    x.pcs <- x.pca$x[,1:pcnum]
    colnames(x.pcs) <- gsub("PC",paste0("PC_",i),colnames(x.pcs))
    x.pcs.l[[i]] <- x.pcs
  }
  ##bind all columns
  x.pcs <- do.call(cbind,x.pcs.l)
  return(x.pcs)
}

### ** function to find correlated genes
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

### * load files
expr <- t(readRDS(expr.file))
expr <- cbind(sample=rownames(expr),data.table(expr))
samples <- expr$sample
expr <- expr[,-c("sample"),with=FALSE]
expr <- as.matrix(expr)
rownames(expr) <- samples
modulegenes <- readLines(module.file)
modulegenes <- gsub("\\..*$","",modulegenes) ## remove extensions of transcript ID
expr.full <- expr
expr <- expr[,intersect(colnames(expr),modulegenes)]

### * PCA
### ** run PCA
expr.pca <- prcomp(expr,retx = TRUE)
### ** find all PCs>1% variance
expr.pca.importance.0 <- data.table(t(summary(expr.pca)$importance))
names(expr.pca.importance.0) <- c("SD","V","PropV")
expr.pca.importance <- expr.pca.importance.0[V>=0.01]
npcs <- nrow(expr.pca.importance)
### ** if no PCs > 1%
if(npcs==0){
  expr.pca.importance <- expr.pca.importance.0[1,]
}
npcs1 <- nrow(expr.pca.importance) ##for random sampling

### * Random genes PCA
### ** generate NULL PCs
glength <- length(modulegenes)
random.pcs <- findpc(n=glength,pcnum=npcs1)
fwrite(expr.pca.importance,paste0(out.file,".importance"),sep="\t",na="NA")
###
### * write PCs
### ** prepare PC file
expr.pcs <- expr.pca$x[,1:npcs]
##align sample order in random pcs
random.pcs <- random.pcs[rownames(expr.pcs),]
expr.pcs <- cbind(expr.pcs,random.pcs)
pc.sample <- rownames(expr.pcs)
expr.pcs <- as.data.table(expr.pcs)
modulename <- basename(module.file)
names(expr.pcs) <- paste0(modulename,"_",names(expr.pcs))

### ** write PC 
if(npcs==0){
  df <- cbind(sample=pc.sample,expr.pcs[,1,with=FALSE])
  fwrite(df,paste0(out.file,0,".expr"), na="NA")
  df.pheno <- cbind(sample=pc.sample,df)
  fwrite(df.pheno,paste0(out.file,0,".pheno"), na="NA", sep="\t", col.names=FALSE)
  siggenes <- findgenes(expr,df)
  writeLines(siggenes,paste0(out.file,0,".sigGenes"))
} else {
  for(i in 1:ncol(expr.pcs)){
    df <- cbind(sample=pc.sample,expr.pcs[,i,with=FALSE])
    pcname <- gsub("^.*PC","",names(expr.pcs)[i])
    fwrite(df,paste0(out.file,pcname,".expr"), na="NA")
    df.pheno <- cbind(sample=pc.sample,df)
    fwrite(df.pheno,paste0(out.file,pcname,".pheno"), na="NA", sep="\t", col.names=FALSE)
    siggenes <- findgenes(expr,df)
    writeLines(siggenes,paste0(out.file,pcname,".sigGenes"))
  }
}


### * org mode specific
### Local Variables:
### eval: (orgstruct-mode 1)
### orgstruct-heading-prefix-regexp: "### "
### End:
