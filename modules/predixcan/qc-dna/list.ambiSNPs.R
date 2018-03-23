#!/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

df.file <- args[1]
output <- args[2]

library(data.table)

df <- fread(df.file)
df$alleles <- paste0(df$V4,df$V5)
df <- df[!alleles %in% c("AT","TA","GC","CG")]

fwrite(df[,"V3",with=FALSE],output,col.names = FALSE,na="NA")
