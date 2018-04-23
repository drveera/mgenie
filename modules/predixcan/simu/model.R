#!/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
gene_name <- args[1]
level_noise <- as.numeric(args[2])

suppressMessages(library(glmnet))
suppressMessages(library(methods))
suppressMessages(library(dplyr))
suppressMessages((library(reshape2)))


gene_annot <- readRDS('/sc/orga/projects/roussp01a/Wen/PredictDB/PredictDBPipeline/data/intermediate/annotations/gene_annotation/gencode.v27.unified.RDS')

                                        # get gene's chromosome and other information

chr_loc=gene_annot[gene_annot$gene_id==gene_name,]
chr_loc=chr_loc[1,]
# load all genotypes on that chromosome
genotype <- read.table(paste0('/sc/orga/projects/roussp01a/Wen/PredictDB/PredictDBPipeline/data/intermediate/genotypes/sim_genotype.chr',chr_loc$chr,'_o.txt'), header = TRUE, row.names = 'rsidID', stringsAsFactors = FALSE)
genotype <- genotype[,-1]
ge1=genotype[,1:2]
genotype <- t(genotype)
gene_annot <- subset(gene_annot, gene_annot$chr == chr_loc$chr)
#snp_annot<-read.table('/sc/orga/projects/roussp01a/Wen/MetaXcan/MetaXcan/prepare/SNP_ann.chr12.txt',header=TRUE,stringsAsFactors=FALSE)
snp_annot <- readRDS(paste0('/sc/orga/projects/roussp01a/Wen/PredictDB/PredictDBPipeline/data/intermediate/annotations/snp_annotation/CMC.snps.anno.chr',chr_loc$chr,'.RDS'))
 rownames(snp_annot)<-snp_annot$rsid
snp_prior<-readRDS(paste0('/sc/orga/projects/roussp01a/Wen/DLPF_annoprior/DLPFCBRN_Anno_prior_rsid.chr',chr_loc$chr,'.RDS'))
rownames(gene_annot) <- gene_annot$gene_id
# weights from maxSNP of eQTL results
eQTL <- read.table(paste0('/sc/orga/projects/epigenAD/coloc/data/CMC/eQTLs/matrixEQTL/DLPFC/chr',chr_loc$chr,'_all_pval.tab'),header=TRUE,stringsAsFactors=FALSE)
geneeQTL <- eQTL[eQTL$ProbeID==gene_name,]

R2_file=matrix(0,100,2)
gene=gene_name
n_k_folds=10
# use the same groupings for both GENET and PrediXcan 
groupid<-read.table('/sc/orga/projects/roussp01a/Wen/Simulation/fixgrouping.txt',header=FALSE,stringsAsFactors=FALSE)
window=1e6
geneinfo <- gene_annot[gene,]
start <- geneinfo$start - window
end <- geneinfo$end + window
cissnps <- subset(snp_annot, snp_annot$pos >= start & snp_annot$pos <= end)
cissnp1index=intersect(colnames(genotype),cissnps$rsid)
rownames(snp_prior)<-snp_prior$RSID_dbSNP137
index1<-intersect(cissnp1index,snp_prior$RSID_dbSNP137)
cisgenos <- genotype[,cissnp1index, drop = FALSE]
prior_penal<-snp_prior[cissnp1index,]
rownames(prior_penal)<-cissnp1index
prior_penal[,3]=rep(0,length(prior_penal[,3]))
prior_penal[index1,3]=snp_prior[index1,3]
cm <- colMeans(cisgenos, na.rm= TRUE)
minorsnps <- subset(colMeans(cisgenos), cm > 0 & cm < 2)
minorsnps <- names(minorsnps)
cisgenos <- cisgenos[,minorsnps, drop = FALSE]
prior_penal<-prior_penal[minorsnps,]
prior_penal<-prior_penal[,3]
pfac=prior_penal
# simulation value
      CMCx2=4.69250216
     # real max prior
     x2=4.75039291
     y2=0.7
     x1=1.9*x2/CMCx2
     y1=0.7
 for (j in 1:length(pfac)) {
          pfac[j]=Re((1+y2-2*y1)/((x2-2*x1)^2)*(2*x1^2+pfac[j]*(x2-2*x1)-2*x1*(x1^2+pfac[j]*(x2-2*x1))^0.5)+(y1-1)/(x2-2*x1)*(-2*x1+2*(x1^2+pfac[j]*(x2-2*x1))^0.5)+1)
          }

#pfac=1/(1+0.33*(pfac)^(0.1))
alpha=0.5

ge1[,2]=rep(0,length(ge1[,2]))

ind<-intersect(geneeQTL$SNPID,colnames(cisgenos))
geneeQTL=geneeQTL[geneeQTL$SNPID %in% ind,]
selectP=min(geneeQTL$PVAL)
selectS=geneeQTL[geneeQTL$PVAL==selectP,]

ge1[selectS$SNPID,2]=selectS$beta
X=genotype

G1exp= X%*%ge1[,2]

for (i in 1:100){
noise = rnorm(n = length(G1exp),mean = 0, sd = 0.3)
G1exp= G1exp+level_noise*noise
exppheno <- scale(G1exp, center = TRUE, scale = TRUE)
# without priors
fit <- tryCatch(cv.glmnet(as.matrix(cisgenos),as.vector(exppheno), nfolds = n_k_folds, alpha = 0.5, type.measure='mse', foldid = groupid$V1, keep = TRUE),
                      error = function(cond) {message('Error'); message(geterrmessage()); list()})
best_lam_ind <- which.min(fit$cvm)
n_fit=fit$nzero[best_lam_ind]
          n_samples=500
          # Discuss with Yungil: use R2 of correlation between predicted and observed expression values for all samples:
          y_all=predict(fit, as.matrix(cisgenos), s = 'lambda.min')
          # correlation R2
          corr_R=ifelse(sd(y_all) != 0, cor(y_all, exppheno),0)
          corr_R2=corr_R**2
          # adjusted correlation R2
          adj_R2=1-(1-corr_R2)*(n_samples-1)/(n_samples-1-n_fit)
          # the second column: ENet
R2_file[i,2] <- adj_R2
# with priors
fit <- tryCatch(cv.glmnet(as.matrix(cisgenos),as.vector(exppheno),penalty.factor=pfac, nfolds = n_k_folds, alpha = 0.5, type.measure='mse', foldid = groupid$V1, keep = TRUE),
                      error = function(cond) {message('Error'); message(geterrmessage()); list()})
best_lam_ind <- which.min(fit$cvm)
n_fit=fit$nzero[best_lam_ind]
          # Discuss with Yungil: use R2 of correlation between predicted and observed expression values for all samples:
          y_all=predict(fit, as.matrix(cisgenos), s = 'lambda.min')
          # correlation R2
          corr_R=ifelse(sd(y_all) != 0, cor(y_all, exppheno),0)
          corr_R2=corr_R**2
          # adjusted correlation R2
          adj_R2=1-(1-corr_R2)*(n_samples-1)/(n_samples-1-n_fit)
# the first column WENet
R2_file[i,1] <- adj_R2
}
print(head(R2_file))
outfile = paste0('/sc/orga/projects/roussp01a/Wen/Simulation/',gene_name,'/adjR2_wenet_enet_',level_noise,'.txt')
write.table(R2_file,outfile,quote=FALSE,col.names=FALSE,row.names=FALSE)

