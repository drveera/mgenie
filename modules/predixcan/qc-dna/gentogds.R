#!/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

genfile <- args[1]
samplefile <- args[2]
outputfile <- args[3]

library(data.table)
library(SNPRelate)
library(GWASTools)

##prepare scadf
sampledf <- fread(samplefile)
sampledf$id <- paste(sampledf$ID_1,sampledf$ID_2)
scandf <- data.frame(sampleID = sampledf$id, scanID = sampledf$ID_2, stringsAsFactors = FALSE)[-1,]

##chromosome
##bim <- fread(genfile,header=FALSE,select=1:2,colClasses = "numeric")
bim <- fread(paste0("cut -d ' ' -f 1,2 ", genfile), header=FALSE)
names(bim) <- c("chr","rsid")
chr <- bim$chr

##convert
imputedDosageFile(input.files=c(genfile,samplefile),
                  filename=outputfile,
                  chromosome = chr,
                  scan.df = scandf)



