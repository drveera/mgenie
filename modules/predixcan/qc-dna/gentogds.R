#!/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

genfile <- args[1]
##samplefile <- args[2]
samplefile <- gsub(".gen$",".sample",genfile)
outputfile <- args[2]

library(data.table)
library(SNPRelate)
library(GWASTools)

##prepare scadf
sampledf <- fread(samplefile)
sampledf$id <- paste(sampledf$ID_1,sampledf$ID_2)
scandf <- data.frame(sampleID = sampledf$id, scanID = sampledf$ID_2, stringsAsFactors = FALSE)[-1,]

##chromosome
chr <- fread(genfile,header=FALSE,select=1,colClasses = "numeric")
chr <- unlist(chr)

##convert
imputedDosageFile(input.files=c(genfile,samplefile), filename=outputfile,chromosome = chr,
                  scan.df = scandf)

