#!/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

gds.file <- args[1]
snpannot.file <- args[2]
genes.file <- args[3]
if(genes.file=="NA"){
    genes.file  <- NA
}
db.file <- args[4]
outfile.gds <- args[5]
outfile.annot <- args[6]

##test arguments
##gds.file <- "merged.gds"
##snpannot.file <- "merged.gds.annot.RDS"
##genes.file <- "sample.gene"
##db.file <- "genie_train/star.aor/star.aor_finalfiles/model.db"
##outfile.gds <- "testout.gds"
##outfile.annot <- "testout.gds.annot.RDS"


library(SNPRelate)
library(GenomicRanges)
library(data.table)
library(dplyr)
library(dbplyr)
library(DBI)
library(RSQLite)
library(GWASTools)

print(gds.file)
print(snpannot.file)
print(genes.file)
print(db.file)
snpannot <- readRDS(snpannot.file)

db <- dbConnect(SQLite(), db.file)
wts <- tbl(db, "weights")
wts <- data.table(as.data.frame(wts))

if(!is.na(genes.file)){
    genes <- fread(genes.file, header=FALSE)
    wts$gene <- gsub("\\..*$","",wts$gene)
    qgenes <- unlist(genes[,1])
    qgenes <- gsub("\\..*$","",qgenes)
    wts <- wts[gene %in% qgenes]
}


snpannot.sub <- snpannot[snpannot$rsid %in% wts$rsid]
saveRDS(snpannot.sub,outfile.annot)
index <- as.data.frame(snpannot.sub)$index

gdsSubset(gds.file, outfile.gds, sample.include = NULL,
          snp.include = index)
