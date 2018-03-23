#!/bin/env Rscript


args <- commandArgs(trailingOnly = TRUE)

bim1.file <- args[1]
bim2.file <- args[2]
outfile <- args[3]

library(data.table)

bim1 <- fread(bim1.file)
bim2 <- fread(bim2.file)

bim1$id <- with(bim1, paste0(V2,V5,V6))
bim2$id <- with(bim2,paste0(V2,V5,V6))
bim3 <- bim1[id %in% bim2$id]

fwrite(bim3[,c("V2"),with=FALSE],outfile,col.names = FALSE,sep="\t",na="NA")
