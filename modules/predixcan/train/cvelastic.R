#!/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
library(glmnet)
library(data.table)
library(GWASTools)
library(SNPRelate)
library(GenomicRanges)
library(dplyr)


##function
cvElastic <- function(gene,geneid,snp,
                      nfolds=10,alpha=0.5){
  set.seed(1)
  groupid <- sample(1:10,length(gene),replace=TRUE)
  fit1 <- cv.glmnet(x=snp,y=gene,nfolds=nfolds,
                    alpha=alpha,foldid = groupid,
                    keep=TRUE,parallel=FALSE)
  fit1.df <- data.table(cvm=fit1$cvm, lambda=fit1$lambda,index=1:length(fit1$cvm))
  best.lambda <- fit1.df[cvm==min(cvm)]
  ##extract betas
  allbetas <- fit1$glmnet.fit$beta[,best.lambda$index]
  allbetas[allbetas==0.0] <- NA
  best.betas <- allbetas[!is.na(allbetas)]
  if(length(best.betas)>2){
      ##names(best.betas) <- row.names(allbetas)[!is.na(allbetas)]
      res <- data.table(rsid=names(best.betas),weight=best.betas)
    res$gene <- geneid
    res$error <- NA
      ##take only useful predictors
      ##allsnps <- data.table(snp=colnames(snp),index=paste0("V",1:ncol(snp)))
      ##usefulsnps <- allsnps[index %in% names(best.betas)]$snp
      snp.new <- snp[,names(best.betas)]

      ##second fit with only the useful snps
      fit2 <- cv.glmnet(x=snp.new, y=gene,
                        nfolds = nfolds, alpha=alpha,
                        foldid=groupid,
                        keep=TRUE,parallel=FALSE)
      fit2.df <- data.table(cvm=fit2$cvm,lambda=fit2$lambda,index=1:length(fit2$cvm))
      best.lambda2 <- fit2.df[cvm==min(cvm)]
      ##global R2
      predicted.gene <- predict(fit2, s=best.lambda2$lambda, newx=snp.new)
      g.rsq <- summary(lm(gene~predicted.gene))$r.squared

      ##true R2
      fold.predicted.gene <- fit2$fit.preval[,best.lambda2$index]
      rsq <- summary(lm(gene~fold.predicted.gene))$r.squared
      res2 <- data.table(gene=geneid,alpha=alpha,
                         n_snps_in_window=ncol(snp),
                         n.snps.in.model=ncol(snp.new),
                         allfolds_R2= rsq,
                         overall_R2 = g.rsq ,
                         error=NA)
      return(list(res,res2))
  } else {
    res2 <- data.table(gene=geneid,alpha=alpha,
                       n_snps_in_window=ncol(snp),
                       n.snps.in.model=length(best.betas),
                       allfolds_R2 = NA,
                       overall_R2 = NA,
                       error=NA)
    res <- data.table(rsid=NA,weight=NA,gene=geneid)
    return(list(res,res2))
  }
}


### codes
gds.file <- args[1]
snpannot.file <- args[2]
genes.file <- args[3]
expr.file <- args[4]
output1 <- args[5]
output2 <- args[6]

####test arguments
##gds.file <- "merged.gds"
##snpannot.file <- "merged.gds.annot.RDS"
##genes.file <- "sample.gene"
##expr.file <- "expression/STARNET.SF.expr.txt.formatted"
##output1 <- "testoutput1"
##output2 <- "testoutput2"

genes <- fread(genes.file, header=FALSE)
names(genes) <- c("gene","chr","start","end")
genes$start <- genes$start - 1000000
genes$end <- genes$end + 1000000

snpannot <- readRDS(snpannot.file)

expr <- fread(expr.file,header=TRUE)
###check later
genes <- genes[gene %in% names(expr)]
#####
expr <- expr[,c("sample",genes$gene),with=FALSE]

genesmodel1 <- list()
genesmodel2 <- list()
for(i in 1:nrow(genes)){
  gene <- genes[i]
  cat("training gene:",gene$gene,"\n")
  gene.ranges <- with(gene, GRanges(seqnames = chr, IRanges(start=start,end=end), gene=gene))
  snpannot.sub <- subsetByOverlaps(snpannot,gene.ranges)
  print(snpannot.sub)
  if(nrow(as.data.frame(snpannot.sub)) < 2){
    r <- data.table(gene=gene$gene,error = "genos is null")
    gene.model <- list(r,r)
  } else {
    index <- as.data.frame(snpannot.sub)$index
    tempgds <- paste0(".temp.",gene$gene,basename(output1),".gds")
    gdsSubset(gds.file,tempgds, snp.include = index)
    mgds <- snpgdsOpen(tempgds)
    genos <- read.gdsn(index.gdsn(mgds, "genotype"))
    snpids <- read.gdsn(index.gdsn(mgds,"snp.rs.id"))
    ##snpids <- as.data.frame(snpannot.sub)$rsid
    rsids <- gsub(",.*$","",snpids)
    chromosomes <- gsub("^.*,","",snpids)
    sampleids <- read.gdsn(index.gdsn(mgds,"sample.id"))
    alleles <- read.gdsn(index.gdsn(mgds,"snp.allele"))
    alleles <- as.data.frame(do.call(rbind,strsplit(alleles, split="/")))
    names(alleles) <- c("ref_allele","eff_allele")
    alleles$rsid <- rsids
    closefn.gds(mgds)
    rownames(genos) <- sampleids
    colnames(genos) <- rsids
    file.remove(tempgds)
    expr1 <- expr[,gene$gene,with=FALSE]
    expr1 <- unlist(expr1)
    names(expr1) <- expr$sample
    genos <- genos[rownames(genos)%in% names(expr1),]
    expr1 <- expr1[rownames(genos)]
    gene.model <- tryCatch(cvElastic(gene=expr1,snp=genos,geneid=gene$gene),
                           error=function(e) {
                             r <- data.table(gene=gene$gene,error=paste0(e))
                             return(list(r,r))
                           } )
    if(!all(is.na(gene.model[[1]]$rsid))) gene.model[[1]] <- merge(gene.model[[1]],alleles,by="rsid")
  }
  genesmodel1[[i]] <- gene.model[[1]]
  genesmodel2[[i]] <- gene.model[[2]]
}

genesmodel1 <- do.call(bind_rows,genesmodel1)
genesmodel2 <- do.call(bind_rows,genesmodel2)

fwrite(genesmodel1,output1,na="NA")
fwrite(genesmodel2,output2,na="NA")
