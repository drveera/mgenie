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
