#!/bin/env Rscript


args <- commandArgs(trailingOnly = TRUE)

gds.file <- args[1]
outfile <- args[2]

library(SNPRelate)
library(gdsfmt)
library(GenomicRanges)

mgds <- snpgdsOpen(gds.file)

rsid <- read.gdsn(index.gdsn(mgds, "snp.rs.id"))
chrom <- gsub("^.*,","",rsid)
rsid <- gsub(",.*$","",rsid)
position <- read.gdsn(index.gdsn(mgds,"snp.position"))
index <- read.gdsn(index.gdsn(mgds,"snp.id"))

gds.ranges <- GRanges(seqnames = chrom, IRanges(start=position,end=position+1),index=index, rsid = rsid)

saveRDS(gds.ranges,outfile)




