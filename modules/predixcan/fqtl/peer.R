#!/bin/env Rscript

###/hpc/packages/minerva-common/R/3.3.1/lib64/R/bin/R

args <- commandArgs(trailingOnly = TRUE)

expr.file <- args[1]
covar.file <- args[2]
outfile1 <- args[3]
outfile2 <- args[4]
pfactors <- as.numeric(args[5])

library(peer)

##pmat <- readRDS("genes.peer.matrix.RDS")
pmat <- readRDS(expr.file)
##covar <- readRDS("covar.peer.RDS")
covar <- readRDS(covar.file)
##covar <- covar[rownames(pmat),]

pmodel <- PEER()
PEER_setPhenoMean(pmodel,pmat)
PEER_setCovariates(pmodel,covar)
dim(PEER_getPhenoMean(pmodel))
PEER_setNk(pmodel,pfactors)
PEER_getNk(pmodel)

PEER_update(pmodel)
##saveRDS(pmodel,"peer.model.RDS")
saveRDS(pmodel,outfile1)

pfactors <- PEER_getX(pmodel)
rownames(pfactors) <- rownames(pmat)
##saveRDS(pfactors,"pfactors.RDS")
saveRDS(pfactors, outfile2)

#weights <- PEER_getW(pmodel)
#prcision <- PEER_getAlpha(pmodel)
#dim(prcision)
#residuals <- PEER_getResiduals(pmodel)



