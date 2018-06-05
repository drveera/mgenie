#!/bin/env Rscript

### * input
args <- commandArgs(trailingOnly = TRUE)
expr.file <- args[1] ## voom object from previous step
covar.file <- paste0(expr.file, ".covar")
outfile1 <- args[2]
outfile2 <- args[3]
pfactors <- as.numeric(args[4])
outfiledge <- args[5]

### * library
library(peer)
library(edgeR)
library(limma)
library(data.table)


### * read files
##expr <- t(readRDS(expr.file)$E)
expr <- readRDS(expr.file)
covar <- fread(covar.file)


### * functions
createDm <- function(x){
  x <- x[,-c("id"),with=FALSE]
  x.classes <- summary.default(x)[,"Mode"]
  if (length(x.classes)==1){
    names(x.classes) <- names(x)
  }
  x.classes.numeric <- x.classes[x.classes=="numeric"]
  x.classes.char <- x.classes[x.classes=="character"]
  if(length(x.classes.char)>0){
    x1 <- x[,names(x.classes.char),with=FALSE]
    x2 <- x[,names(x.classes.numeric),with=FALSE]
    fm <- as.formula(paste0("~",paste0(names(x1),collapse = "+")))
    x1.dm <- model.matrix(fm,data=x1)
    print(dim(x1.dm))
    print(dim(x))
    return(cbind(x1.dm,as.matrix(x2)))
  } else {
    x2 <- x[,names(x.classes.numeric),with=FALSE]
    return(as.matrix(x2))
  }
}


### * process data 
### ** subset samples common to expr and covar
covar <- covar[complete.cases(covar)]
ids.merge <- intersect(rownames(expr),covar$id)
if(length(ids.merge)==0) stop("no samples from covar file match with rownames of expression file")
covar <- covar[id %in% ids.merge]
expr <- expr[covar$id,] ## this way the order of samples in expr is same as covar

### ** run dge, save dge for next step and run voom for peer
exprt <- t(expr)
dge <- DGEList(exprt)
cutoff <- 0.2
genestokeep <- rowSums(cpm(dge)>1) >= (ncol(dge)*cutoff)
dge <- dge[genestokeep,]
dge <- calcNormFactors(dge, method="TMM")
saveRDS(dge,outfiledge)
expr.voom <- voom(dge,plot=FALSE)
expr <- t(expr.voom$E)
print(table(rownames(expr) == covar$id))

### ** convert covar to matrix
print(head(covar))
covar.mat <- createDm(covar)

### * PEER
pmodel <- PEER()
PEER_setPhenoMean(pmodel,expr)
PEER_setCovariates(pmodel,covar.mat)
dim(PEER_getPhenoMean(pmodel))
PEER_setNmax_iterations(pmodel, 1000)
PEER_setNk(pmodel,pfactors)
PEER_getNk(pmodel)
PEER_update(pmodel)
### * output 
saveRDS(pmodel,outfile1)
pfactors <- PEER_getX(pmodel)
print(head(pfactors))
rownames(pfactors) <- rownames(expr)
##saveRDS(pfactors,"pfactors.RDS")
saveRDS(pfactors, outfile2)

### * org mode specific
### Local Variables:
### eval: (orgstruct-mode 1)
### orgstruct-heading-prefix-regexp: "### "
### End:


