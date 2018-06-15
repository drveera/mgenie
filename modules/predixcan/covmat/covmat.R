#!/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

gds.file <- args[1]
snpannot.file <- args[2]
genes.file <- args[3]
db.file <- args[4]
outfile <- args[5]

##gds.file <- "1kgsub.gds"
##snpannot.file <- "1kgsub.gds.annot.RDS"
##genes.file <- "cmcgene000"
##db.file <- "gtex_v7_Adipose_Subcutaneous_imputed_europeans_tw_0.5_signif.db"

if(FALSE){
  gds.file <- "1kg.gds"
  snpannot.file <- "1kg.gds.annot.RDS"
  genes.file <- "/hpc/users/xrajagv01/va-biobank/Veera/analysis/modules/GTEx/modulegenes/Adipose...Visceral._ob_Omentum_cb_.GTEX.expression.RDS.0.pcovar.residuals_color_beige.module.genes"
  db.file <- "module.db"
  outfile <- "testoutput"
}

library(SNPRelate)
library(GenomicRanges)
library(data.table)
library(dplyr)
library(dbplyr)
library(DBI)
library(RSQLite)
library(GWASTools)
library(reshape)





db <- dbConnect(SQLite(), db.file)
wts <- tbl(db, "weights")
wts <- data.table(as.data.frame(wts))
##wts$gene1 <- gsub("\\..*$","",wts$gene)
##print(head(genes))
##print(head(wts$gene))
if (!genes.file=="NA"){
  genes <- fread(genes.file, header=FALSE)
  genes <- unlist(genes[,1])
  ##genes <- gsub("\\..*$","",genes)
  genes <- intersect(genes,wts$gene)
} else {
  genes <- wts$gene
}


snpannot <- readRDS(snpannot.file)
##mgds <- snpgdsOpen(gds.file)

qgene <- genes[3]
cvmat <- function(qgene,wts,snpannot,gds.file){
  wts.sub <- wts[gene==qgene]
  fgene <- wts.sub$gene[1] ## gene name to write to result
  qrsids <- wts.sub$rsid
  snpannot.sub <- snpannot[snpannot$rsid %in% qrsids,]
  index.sub <- as.data.frame(snpannot.sub)$index
  if(length(index.sub)<2){
    return(NULL)
  }
  tempgds <- paste0(".temp.",qgene,".gds")
  gdsSubset(gds.file,tempgds, snp.include = index.sub)
  mgds <- snpgdsOpen(tempgds)
  genos <- read.gdsn(index.gdsn(mgds,"genotype"))
  snps <- read.gdsn(index.gdsn(mgds,"snp.rs.id"))
  closefn.gds(mgds)
  file.remove(tempgds)
  snps <- gsub(",.*$","",snps)
##  print(head(genos))
##  print(class(genos))
  colnames(genos) <- snps
  print(dim(genos))
  covmat <- cov(genos)
  covmat.melt <- melt(covmat)
  covmat.melt <- data.table(covmat.melt)
  names(covmat.melt) <- c("RSID1","RSID2","VALUE")
  covmat.melt$GENE <- fgene
  covmat.melt <- covmat.melt[,c("GENE","RSID1","RSID2","VALUE")]
  return(covmat.melt)
}

results <- list()
for(i in genes){
  cat("working on gene:",i,"\n")
  results[[i]] <- cvmat(i,wts,snpannot,gds.file)
}
covmat.merged <- do.call(rbind,results)
##print(head(covmat.merged))
##print(class(covmat.merged))
if(is.null(covmat.merged)){
  system(paste0("touch ",outfile))
} else {
  fwrite(covmat.merged,outfile,sep="\t",na="NA")
}
