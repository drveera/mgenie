#!/bin/env Rscript


args <- commandArgs(trailingOnly = TRUE)

sex.file <- args[1]
outfile <- args[2]

library(data.table)

sex <- fread(sex.file)

sex <- sex[STATUS!="OK"]

fwrite(sex,outfile,sep="\t",na="NA")







