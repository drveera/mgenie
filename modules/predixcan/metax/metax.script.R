#!/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

gwas.file <- args[1]
db.file <- args[2]
print(db.file)
cov.file <- gsub("db$","cov.gz",db.file)
genes.file <- args[3]
genes <- readLines(genes.file)
print(genes.file)
out.file <- args[4]

library(metaxcanr)
library(data.table)
res <- metaxcan(gwas.file = gwas.file,
                db.file = db.file, snpcov.file = paste0("zcat ",cov.file),
                genes = genes)

if(is.null(res)){
  file.create(out.file)
}else {
  fwrite(res,out.file,sep="\t",na="NA")
}
