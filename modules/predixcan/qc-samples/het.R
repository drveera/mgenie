#!/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

het.file <- args[1]
fam.file <- args[2]
outfile <- args[3]

library(data.table)

het <- fread(het.file)
het <- het[F>0.2]

fam <- fread(paste0(fam.file,".fam"))
fam <- fam[V2 %in% het$IID]

fwrite(fam,outfile,sep="\t",na="NA")

