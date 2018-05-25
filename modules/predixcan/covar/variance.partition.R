#!/bin/env Rscript

### * arguments
args <- commandArgs(trailingOnly = TRUE)
expr.file <- args[1]
covar.file <- args[2]
out.file <- args[3]

### * libraries
library(data.table)
library(variancePartition)
library(limma)
library(doParallel)
library(parallel)

### * read files
expr <- readRDS(expr.file)
expr_E <- expr$E
dim(expr)
##expr <- t(expr)
covar <- fread(covar.file)
covar <- covar[id %in% colnames(expr_E)]
rownames(covar) <- covar$id
dim(covar)
##dim(expr)


### * varianc partition

### ** create formula
covar <- covar[,-"id",with=FALSE]
covar.names <- names(covar)
### column classes
cnames.class <- summary.default(covar)[,"Mode"]
##get variable names that are numeric
cnames.class.numeric <- names(cnames.class[cnames.class=="numeric"])
##get variable names that are character
cnames.class.chr <- names(cnames.class[cnames.class=="character"])
if(length(cnames.class.chr)>0){
  cnames.class.chr <- paste0("(1|",cnames.class.chr,")")
} else {
  cnames.class.chr <- ""
}
##make formula
fvars <- c(cnames.class.chr,cnames.class.numeric)
fmula <- as.formula(paste0("~",paste0(fvars,collapse = "+")))

### ** set up cluster
cl <- makeCluster(8)
registerDoParallel(cl)
print(fmula)
print(dim(covar))
print(table(covar$bsscode))
### ** variance parition
vpart <- fitExtractVarPartModel(expr,fmula,covar)
saveRDS(vpart,out.file)
vpart.plot <- plotVarPart(vpart)
vpart.plot1 <- vpart.plot + theme(axis.text.x=element_text(angle=90, hjust=1, size=7))
ggsave(plot=vpart.plot1, filename=paste0(out.file,"plot.pdf"))
##ggsave("~/bridge/vpart.pdf")

### * org mode specific
### Local Variables:
### eval: (orgstruct-mode 1)
### orgstruct-heading-prefix-regexp: "### "
### End:
