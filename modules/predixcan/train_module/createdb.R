#!/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

wts.file <- args[1]
extra.file <- args[2]
db.out <- args[3]


library(dplyr)
library(dbplyr)
library(DBI)
library(RSQLite)
library(data.table)


wts <- fread(wts.file)
wts$weight <- as.numeric(wts$weight)
wts <- wts[!is.na(weight)]
wts <- wts[,c("gene","rsid","weight","ref_allele","eff_allele"),with=FALSE]

extra <- fread(extra.file)
extra$alpha <- as.numeric(extra$alpha)


mydb <- dbConnect(SQLite(),db.out,create=TRUE)

copy_to(mydb,df=wts,name="weights",temporary = FALSE)
copy_to(mydb,df=extra,name="extra",temporary = FALSE)

chunk <- 50
n <- nrow(extra)
r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
extra.chunks <- split(extra,r)

pb <- txtProgressBar(min = 0, max = length(extra.chunks), style = 3)
for(i in 1:length(extra.chunks)){
    extra.sub  <- extra.chunks[[i]]
    wts.sub  <- wts[gene %in% extra.sub$gene_id]
    db.out.sub  <- paste0(db.out,".chunk",i,".db")
    mydb <- dbConnect(SQLite(),db.out.sub,create=TRUE)
    copy_to(mydb,df=wts.sub,name="weights",temporary = FALSE)
    copy_to(mydb,df=extra.sub,name="extra",temporary = FALSE)
    setTxtProgressBar(pb, i)
}
close(pb)



