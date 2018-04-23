#!/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

flist <- args[1]
out.file <- args[2]

library(data.table)
library(dplyr)
flist <- readLines(flist)
df <- list()
for(i in flist) df[[i]] <- fread(i)
df <- do.call(bind_rows,df)

fwrite(df,out.file,na="NA")
