#!/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

gwas.file <- args[1]
db.file <- args[2]
cov.file <- paste0(db.file,".covmat")
out.file <- args[3]
genes.file <- args[4]
print(genes.file)
if(!is.na(genes.file)){
  genes <- readLines(genes.file)
} else {
  genes <- NA
}

library(metaxcanr)
library(data.table)
res <- metaxcan(gwas.file = gwas.file,
                db.file = db.file, snpcov.file = cov.file,
                genes = genes, ncores=4)

if(is.null(res)){
  file.create(out.file)
}else {
  fwrite(res,out.file,sep="\t",na="NA")
}
