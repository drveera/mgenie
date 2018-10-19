#!/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

flist <- args[1]
out.file <- args[2]

library(data.table)
library(dplyr)
flist <- readLines(flist)
df <- list()
fread_special <- function(x){
  ##xdf <- fread(x)
  xdf <- read.table(x,header = TRUE, stringsAsFactors = FALSE, sep=",")
  xdf <- data.table(xdf)
  if("ref_allele" %in% names(xdf)){
    xdf <- fread(x,colClasses = list(character=5:6))
  }
  return(xdf)
}
for(i in flist) df[[i]] <- fread_special(i)
df <- do.call(bind_rows,df)

fwrite(df,out.file,na="NA",sep="\t")
