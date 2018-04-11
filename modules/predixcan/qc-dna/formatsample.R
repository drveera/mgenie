#!/bin/env Rscript


args <- commandArgs(trailingOnly = TRUE)

sample.file <- args[1]
outfile <- args[2]

library(data.table)

sample <- fread(sample.file)

sample$ID_1 <- sample$ID_2

fwrite(sample,outfile,sep="\t",na="NA")


