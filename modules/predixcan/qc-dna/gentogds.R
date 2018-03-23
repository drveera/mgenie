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
print(head(sampledf))
##sampledf$id <- paste(sampledf$ID_1,sampledf$ID_2)
sampledf$id <- sampledf$ID_2
scandf <- data.frame(sampleID=1:nrow(sampledf),scanID=sampledf$id)[-1,]
scandf$scanID <- as.character(scandf$scanID)


##chromosome
chr <- fread(genfile,header=FALSE,select=1,colClasses = "numeric")
chr <- unlist(chr)

##convert
imputedDosageFile(input.files=c(genfile,samplefile), filename=outputfile,chromosome = chr,
                  scan.df = scandf)
