#!/bin/env Rscript


args <- commandArgs(trailingOnly = TRUE)

gds.file <- args[1]
genfile <- args[2]
outfile <- args[3]


library(SNPRelate)
library(gdsfmt)
library(GenomicRanges)
library(data.table)

mgds <- snpgdsOpen(gds.file)

position <- read.gdsn(index.gdsn(mgds,"snp.position"))
alleles <- read.gdsn(index.gdsn(mgds,"snp.allele"))
index <- read.gdsn(index.gdsn(mgds,"snp.id"))

##gen file
##gen <- fread(genfile,select=1:5)
gen <- fread(paste0("cut -d ' ' -f 1-5 ",genfile), header=FALSE)
gen$alleles <- paste0(gen$V4,"/",gen$V5)

##checks
if(all(position == gen$V3)){
  message("all positions match between gen and gds")
} else {
  stop("position mismatch between gen and gds")
}

if(all(alleles == gen$alleles)) {
  message("all alleles match between gen and gds")
} else {
  stop("alleles mismatch between gen and gds")
}

gds.ranges <- GRanges(seqnames = gen$V1,
                      IRanges(start=position,end=position+1),
                      index=index,
                      rsid = gen$V2,
                      A1=gen$V4,
                      A2=gen$V5)

saveRDS(gds.ranges,outfile)





