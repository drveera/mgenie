#!/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

bfile <- args[1]
reffam <- args[2]
outfile <- args[3]
orginalfam <- args[4]

library(SNPRelate)
library(data.table)
library(ggplot2)

##conver to gds
snpgdsBED2GDS(bed.fn = paste0(bfile,".bed"),
              fam.fn = paste0(bfile,".fam"),
              bim.fn = paste0(bfile,".bim"),
              out.gdsfn = paste0(outfile,".gds"))

##open the gds
genofile <- snpgdsOpen(paste0(outfile,".gds"))
##run pca
pca <- snpgdsPCA(genofile,num.thread = 4)
pca.eigen <- pca$eigenvect
pca.eigen <- data.table(pca.eigen)
IID <- pca$sample.id
pca.eigen <- cbind(IID,pca.eigen)

##read referrence samples
fam <- fread(reffam)
refsamples = fam$V1
##read bfile samples
bfam <- fread(paste0(bfile,".fam"))
bfamsamples <- bfam$V1

pca.eigen$popgroup <- ifelse(pca.eigen$IID %in% refsamples,"Referrence","Cases")

##plot1
##p1 <- ggplot(pca.eigen, aes(V1,V2, colour=popgroup)) + geom_point()
##ggsave(plot1)


##subset the dfm
pdfm = pca.eigen[IID %in% refsamples,]
##formula
pca.eigen$ell2sd6 <- (
  ((pca.eigen$V1 - mean(pdfm$V1, na.rm=TRUE))/(6*sd(pdfm$V1,na.rm=TRUE)))^2 +
  ((pca.eigen$V2 - mean(pdfm$V2, na.rm=TRUE))/(6*sd(pdfm$V2, na.rm=TRUE)))^2
) <=1
##pca.eigen$ell2sd6 = pca.eigen$ell2sd6*1

##pca.eigen$popgroup <- with(pca.eigen, ifelse(ell2sd6 & (!IID %in% refsamples),"keep",popgroup))
##pca.eigen$popgroup <- with(pca.eigen, ifelse((!ell2sd6) & (!IID %in% refsamples),"remove",popgroup))

pca.eigen$popgroup <- with(pca.eigen, ifelse(ell2sd6 & (IID %in% bfamsamples),"keep",popgroup))
pca.eigen$popgroup <- with(pca.eigen, ifelse((!ell2sd6) & (IID %in% bfamsamples),"remove",popgroup))

p2 <- ggplot(pca.eigen, aes(V1,V2, colour=popgroup)) + geom_point()
plot2 <- paste0(outfile,".PCA.pdf")
ggsave(plot2)

fwrite(pca.eigen,paste0(outfile,".allsamples.mds"),sep="\t",na="NA")
##pca.eigen1 <- pca.eigen[ell2sd6==TRUE & popgroup!="Referrence",]
pca.eigen1 <- pca.eigen[popgroup="keep",]

##repeat PCA
pca2 <- snpgdsPCA(genofile,num.thread = 4,sample.id=pca.eigen1$IID)
pca2.eigen <- pca2$eigenvect
pca2.eigen <- data.table(pca2.eigen)
IID <- pca2$sample.id
pca2.eigen <- cbind(IID,pca2.eigen)

fwrite(pca2.eigen,paste0(outfile,".mds"),sep="\t",na="NA")
##selected.fam <- data.table(IID=pca2.eigen$IID,FID=pca2.eigen$IID,PHENO=1)
selected.fam <- fread(orginalfam)[V2 %in% pca2.eigen$IID]
fwrite(selected.fam,paste0(outfile,".samples"),sep="\t",col.names = TRUE, na="NA",quote=FALSE)
