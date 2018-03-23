#!/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

gds.file <- args[1]
snpannot.file <- args[2]
genes.file <- args[3]
db.file <- args[4]
output1 <- args[5]
output2 <- args[6]


if(FALSE){
gds.file <- "chrom22.gds"
snpannot.file <- "snp.annot.RDS"
genes.file <- "temp.genes"
db.file <- "ebv.db"
output1 <- "out1"
output2 <- "out2"
}

library(DBI, lib.loc = "~/va-biobank/Veera/Rlibraries/")
library(GWASTools)
library(data.table)
library(SNPRelate)
library(GenomicRanges)
library(dbplyr)
library(dplyr)
library(RSQLite)

snpannot <- readRDS(snpannot.file)
genes <- readLines(genes.file)
db <- dbConnect(SQLite(), db.file)

totalexpr <- list()
for(i in 1:length(genes)){
  ensid <- genes[i]
  ##get weights
  dbwts <- tbl(db,"weights") %>% filter(gene == ensid) %>% collect()
  dbwts <- as.data.frame(dbwts)
  rownames(dbwts) <- dbwts$rsid
  snps.df <- snpannot[snpannot$rsid %in% dbwts$rsid]
  ##this snps.df must have columns rsid,index,A1,A2
  if(nrow(as.data.frame(snps.df))==0){
    mgds <- snpgdsOpen(gds.file)
    sampleids <- read.gdsn(index.gdsn(mgds,"sample.id"))
    closefn.gds(mgds)
    predictedexpr <- matrix(NA,nrow=1,ncol=2504)
    colnames(predictedexpr) <- sampleids
    totalexpr[[i]] <- predictedexpr
  } else {
    ##update dbwts order 
    dbwts <- dbwts[snps.df$rsid,]
    snpwts <- dbwts$weight
    names(snpwts) <- dbwts$rsid
    ##get index
    index <- snps.df$index
    ##subset gds first
    tempgds <- paste0(dirname(output1),"/.temp",ensid,".gds")
    gdsSubset(gds.file,tempgds,snp.include = index)
    mgds <- snpgdsOpen(tempgds)
    genos <- read.gdsn(index.gdsn(mgds,"genotype"))
    snporder <- read.gdsn(index.gdsn(mgds,"snp.id"))
    sampleids <- read.gdsn(index.gdsn(mgds,"sample.id"))
    alleles <- read.gdsn(index.gdsn(mgds,"snp.allele"))
    alleles <- do.call(rbind, strsplit(alleles,split="/"))
    alleles <- data.table(alleles)
    names(alleles) <- c("A1","A2")
    ##genos matrix file rows samples and columns snps
    ##need to  invert the dose and change according to alleles
    ##this is always true, since gds format by default stores the count of first allele (A1)
    ##by subtracting from 2, the dose will represent A2
    ##so, check if A2 == eff_allele 
    genos <- 2 - genos
    ##first check the samples are in columns and genes in rows
    if(length(snpwts) != nrow(genos)){
      genos <- t(genos)
    }
    ##check alleles
    for(i in nrow(genos)){
      if(alleles$A2[i] != dbwts$eff_allele[i]
         & alleles$A2[i] == dbwts$ref_allele[i]
         ){
        genos[i,] <- 2 - genos[i,]
      }
    }
    predictedexpr <- snpwts %*% genos
    colnames(predictedexpr) <- sampleids
    totalexpr[[i]] <- predictedexpr
    ##finally
    ##close and delete the gds
    closefn.gds(mgds)
    file.remove(tempgds)
  }
}


totalexpr1 <- do.call(rbind,totalexpr)
totalexpr <- data.table(totalexpr1)
totalexpr <- cbind(genes,totalexpr)


fwrite(totalexpr,output1,na="NA")
###will change later for now just copy the output1 to output2
fwrite(totalexpr,output2,na="NA")


