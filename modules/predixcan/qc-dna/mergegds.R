#!/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

gds.files.file <- args[1]
out.file <- args[2]


library(gdsfmt)

##gds.files <- list.files(pattern = "*copy.gds$")
gds.files <- readLines(gds.files.file)
gds.list <- list()
for(i in  gds.files) gds.list[[i]] <- openfn.gds(i)

## check if all samples are same
sampleids <- lapply(gds.list, function(x) return(read.gdsn(index.gdsn(x,"sample.id"))))
refids <- sampleids[[1]]
if (!all(sapply(sampleids, function(x) return(all(x==refids))))) stop("sample are not same across gds files")

##snpids <- lapply(gds.list, function(x) return(read.gdsn(index.gdsn(x,"snp.id"))))
rsids <- do.call(c,lapply(gds.list, function(x) return(read.gdsn(index.gdsn(x,"rsid")))))
## snp positions
snp.positions <- do.call(c,lapply(gds.list, function(x) return(read.gdsn(index.gdsn(x,"snp.position")))))
## snp.chromosome
snp.chromosomes <- do.call(c,lapply(gds.list, function(x) return(read.gdsn(index.gdsn(x,"snp.chromosome")))))
## snp alleles
snp.alleles <- do.call(c,lapply(gds.list, function(x) return(read.gdsn(index.gdsn(x,"snp.allele")))))

## check if the all have same names
clist <- list(rsids,snp.positions,snp.chromosomes,snp.alleles)
refnames <- names(clist[[1]])
if(!all(sapply(clist,function(x) all(names(x)==refnames)))) stop("names don't match between elements")
## genotypes
genotypes <- do.call(rbind,lapply(gds.list, function(x) return(read.gdsn(index.gdsn(x,"genotype")))))

## create gds file
f <- createfn.gds(out.file)
add.gdsn(f, "snp.id",1:length(rsids), compress = "LZMA_RA.fast")
add.gdsn(f, "sample.id",refids, compress = "LZMA_RA.fast")
add.gdsn(f, "genotype",genotypes, compress = "LZMA_RA.fast")
add.gdsn(f, "snp.chromosome", snp.chromosomes, compress = "LZMA_RA.fast")
add.gdsn(f, "snp.position", snp.positions, compress = "LZMA_RA.fast")
add.gdsn(f, "snp.allele", snp.alleles, compress = "LZMA_RA.fast")
add.gdsn(f, "rsid",rsids,compress = "LZMA_RA.fast")
closefn.gds(f)

