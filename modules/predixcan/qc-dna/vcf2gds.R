#!/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

vcf.file <-  args[1]
outfile <- args[2]
print(outfile)

library(SNPRelate)

snpgdsVCF2GDS(vcf.file,outfile)

