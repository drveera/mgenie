
args <- commandArgs(trailingOnly = TRUE)

plinkbase <- args[1]
outfile <- args[2]

library(SNPRelate)

bedfile <- paste0(plinkbase,".bed")
bimfile <- paste0(plinkbase,".bim")
famfile <- paste0(plinkbase,".fam")
snpgdsBED2GDS(bedfile, famfile, bimfile, outfile)
