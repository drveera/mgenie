#!/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

gds.file <- args[1]
snpannot.file <- args[2]
genes.file <- args[3]
db.file <- args[4]
output1 <- args[5]
output2 <- args[6]


if(FALSE){
gds.file <- "temp.del.gds"
snpannot.file <- "temp.gds.annot"
genes.file <- "genebatch/genes000"
db.file <- "cmcphase1.nopriors.nogroups.db"
output1 <- "temp.del.out1"
output2 <- "temp.del.out2"
}

##library(DBI, lib.loc = "~/va-biobank/Veera/Rlibraries/")
library(DBI)
library(GWASTools)
library(data.table)
library(SNPRelate)
library(GenomicRanges)
library(dbplyr)
library(dplyr)
library(RSQLite)

snpannot <- readRDS(snpannot.file)
genes <- unlist(fread(genes.file)[,1,with=FALSE])
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
    predictedexpr <- matrix(NA,nrow=1,ncol=length(sampleids))
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
    genos <- as.matrix(read.gdsn(index.gdsn(mgds,"genotype")))
    ##as matrix because if only one SNP, then it loses matrix class
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
    print(length(snpwts))
    print(dim(genos))
    if(length(snpwts) != nrow(genos)){
      genos <- t(genos)
    }
    ##check alleles
    for(j in nrow(genos)){
      if(alleles$A2[j] != dbwts$eff_allele[j]
         & alleles$A2[j] == dbwts$ref_allele[j]
         ){
        genos[j,] <- 2 - genos[j,]
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
totalexpr <- t(totalexpr)
##colnames(totalexpr) <- c("id",genes)
##totalexpr <- cbind(genes,totalexpr)
##genes in columns and samples in rows
totalexpr <- data.table(totalexpr)

fwrite(totalexpr,output1,na="NA",quote=FALSE)

###will change later for now just copy the output1 to output2
fwrite(totalexpr,output2,na="NA",quote=FALSE)


