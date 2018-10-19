#!/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

gds.file <- args[1]
snpannot.file <- args[2]
genes.file <- args[3]
outgds <- args[4]
outannot <- args[5]


if (FALSE){
  gds.file <- "philipshaw.gds"
  snpannot.file <- "philipshaw.gds.annot.RDS"
  genes.file <- "genebatch/genes000"
  outgds <- "temp.del.gds"
  outannot <- "temp.gds.annot"
}

library(GWASTools)
library(SNPRelate)
library(GenomicRanges)
library(data.table)

snpannot <- readRDS(snpannot.file)
genes <- fread(genes.file, header=FALSE)
names(genes) <- c("gene","chr","start","end")
genes$start <- genes$start - 1000000
genes$end <- genes$end + 1000000
genes.ranges <- with(genes, GRanges(seqnames = chr, IRanges(start=start,end=end), gene=gene))
snpannot.sub <- subsetByOverlaps(snpannot,genes.ranges)
saveRDS(snpannot.sub,outannot)
##print(head(snpannot.sub))
index <- as.data.frame(snpannot.sub)$index

##subset gds
gdsSubset(gds.file, outgds, sample.include = NULL,
          snp.include = index)
snpgdsSummary(outgds)

