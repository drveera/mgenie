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
extra <- fread(extra.file)

mydb <- dbConnect(SQLite(),db.out,create=TRUE)

copy_to(mydb,df=wts,name="weights",temporary = FALSE)
copy_to(mydb,df=extra,name="extra",temporary = FALSE)
