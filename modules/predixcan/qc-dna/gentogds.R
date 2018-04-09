#!/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

genfile <- args[1]
##samplefile <- args[2]
samplefile <- gsub(".gen$",".sample",genfile)
outputfile <- args[2]
gdscopy.file <- args[3]

library(data.table)
library(SNPRelate)
library(GWASTools)

##prepare scadf
sampledf <- fread(samplefile)
sampledf$id <- paste(sampledf$ID_1,sampledf$ID_2)
scandf <- data.frame(sampleID = sampledf$id, scanID = sampledf$ID_2, stringsAsFactors = FALSE)[-1,]

##chromosome
bim <- fread(genfile,header=FALSE,select=1:2,colClasses = "numeric")
names(bim) <- c("chr","rsid")
chr <- bim$chr

##convert
imputedDosageFile(input.files=c(genfile,samplefile), filename=outputfile,chromosome = chr,
                  scan.df = scandf)

##temp
mgds <- openfn.gds(outputfile)
f <- createfn.gds(gdscopy.file)
add.gdsn(f, "snp.id",read.gdsn(index.gdsn(mgds,"snp.id")),compress = "LZMA_RA.fast")
add.gdsn(f, "sample.id",read.gdsn(index.gdsn(mgds,"sample.id")),compress = "LZMA_RA.fast")
add.gdsn(f, "genotype",read.gdsn(index.gdsn(mgds,"genotype")),compress = "LZMA_RA.fast", storage="float32")
add.gdsn(f, "snp.chromosome",read.gdsn(index.gdsn(mgds,"snp.chromosome")),compress = "LZMA_RA.fast")
add.gdsn(f, "snp.position",read.gdsn(index.gdsn(mgds,"snp.position")),compress = "LZMA_RA.fast")
add.gdsn(f, "snp.allele",read.gdsn(index.gdsn(mgds,"snp.allele")),compress = "LZMA_RA.fast")
add.gdsn(f,"rsid",bim$rsid,compress = "LZMA_RA.fast")
closefn.gds(f)


