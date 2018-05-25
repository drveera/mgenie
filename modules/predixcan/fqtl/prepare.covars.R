#!/bin/env Rscript

### * arguments
args <- commandArgs(trailingOnly = TRUE)
peerfactors.file <- args[1]
covarout <- args[2]
npeer <- as.numeric(args[3])

### * library
library(data.table)

### * read files
##covar <- fread(covar.file)
peerfactors <- readRDS(peerfactors.file)

### * prepare peer added covar file
peerfactors <- as.data.frame(peerfactors)
peerfactors$id <- rownames(peerfactors)
peerfactors <- peerfactors[,c(ncol(peerfactors),1:(ncol(peerfactors)-1))]
print(head(peerfactors))
pstart <- ncol(peerfactors) - npeer
j=0
for(i in pstart:ncol(peerfactors)){
  df <- peerfactors[,1:i]
  fwrite(df, paste0(covarout,".",j,".pcovar"), sep="\t", na="NA")
  j=j+1
}

### * org mode sepcific
### Local Variables:
### eval: (orgstruct-mode 1)
### orgstruct-heading-prefix-regexp: "### "
### End:
