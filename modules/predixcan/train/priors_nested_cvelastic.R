#!/bin/env Rscript

### * libraries
library(glmnet)
library(data.table)
library(GWASTools)
library(SNPRelate)
library(GenomicRanges)
library(dplyr)

### * arguments
args <- commandArgs(trailingOnly = TRUE)
gds.file <- args[1]
snpannot.file <- args[2]
genes.file <- args[3]
expr.file <- args[4]
output1 <- args[5]
output2 <- args[6]
priordf.file <- args[7]
grouping.file <- args[8]



if(FALSE){
####test arguments
gds.file <- "merged.gds"
snpannot.file <- "merged.gds.annot.RDS"
genes.file <- "test.gene"
##genes.file <- "sample.genes.bed"
##expr.file <- "expression/STARNET.SF.expr.txt.formatted"
##expr.file <- "alvaro.gtex.brain.cerebellum.residuals"
expr.file <- "CMC.expr.txt.formatted"
output1 <- "priortestoutput1"
output2 <- "priortestoutput2"
priordf.file <- "../../../priors/dlpfc/full.anno.priors.txt.formatted"
grouping.file <- "../../../grouping/genoCMCgroup.txt.formatted"
}

### * functions
### ** cvElastic
cvElastic <- function(gene,
                      geneid,
                      snp,pfac,
                      grouping.file,
                      nfolds=10,
                      alpha=0.5){
  set.seed(1)
  if(is.na(grouping.file)){
    groupid <- sample(1:10,length(gene),replace=TRUE)
  } else {
    groupid <- manual_grouping(grouping.file,"grouping",gene)
  }
  fit1 <- cv.glmnet(x=snp,y=gene,
                    nfolds=nfolds,
                    alpha=alpha,
                    foldid = groupid,
                    keep=TRUE,
                    parallel=FALSE,
                    penalty.factor=pfac,
                    type.measure = "mse")
  fit1.df <- data.table(cvm=fit1$cvm, lambda=fit1$lambda,index=1:length(fit1$cvm))
  best.lambda <- fit1.df[cvm==min(cvm)]
  ##extract betas
  allbetas <- fit1$glmnet.fit$beta[,best.lambda$index]
  allbetas[allbetas==0.0] <- NA
  best.betas <- allbetas[!is.na(allbetas)]
  if(length(best.betas)>=2){
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
                        keep=TRUE,
                        parallel=FALSE,
                        penalty.factor=pfac)
      fit2.df <- data.table(cvm=fit2$cvm,lambda=fit2$lambda,index=1:length(fit2$cvm))
    best.lambda2 <- fit2.df[cvm==min(cvm)]
    ##adjusted R2 #@@@2
    a.pred.values <- predict(fit1, s=best.lambda$lambda,newx=snp)
    a.corr.R2 <- ifelse(sd(a.pred.values) != 0, cor(a.pred.values,gene),0)
    a.corr.R2sq <- a.corr.R2**2
    ##a.adj.R2sq <- 1 - (1-corr.R2)*(n_)
    #@@@@

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
    ##nested.cv performance evaluation
    prf.stats <- nested_cv_elastic_net_perf(x=snp,
                                            y=gene,
                                            n_samples = length(gene),
                                            n_train_test_folds = 5,
                                            n_k_folds=10,
                                            alpha=0.5,
                                            pfac=pfac,
                                            grouping.file=grouping.file)
    prf.stats <- do.call(cbind,prf.stats)
    prf.stats <- data.table(prf.stats)
    names(prf.stats) <- c("rmse_avg",
                          "test_R2_avg",
                          "test_R2_sd",
                          "nested_cv_fisher_pval",
                          "rho_avg",
                          "rho_se",
                          "rho_zscore",
                          "pred.perf.R2",
                          "pred.perf.pval")
    res2 <- cbind(res2,prf.stats)
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

### ** small functions
### *** generate_fold_ids
generate_fold_ids <- function(n_samples, n_folds=10) {
  n <- ceiling(n_samples / n_folds)
  fold_ids <- rep(1:n_folds, n)
  sample(fold_ids[1:n_samples])
}

### *** adjust_for_covariates
adjust_for_covariates <- function(expression_vec, cov_df) {
  combined_df <- cbind(expression_vec, cov_df)
  expr_resid <- summary(lm(expression_vec ~ ., data=combined_df))$residuals
  expr_resid <- scale(expr_resid, center = TRUE, scale = TRUE)
  expr_resid
}

### *** calc_R2
calc_R2 <- function(y, y_pred) {
  ##tss <- sum(y**2)
  tss <- sum((y-mean(y))**2)
  rss <- sum((y - y_pred)**2)
  1 - rss/tss
}

### *** calc_corr
calc_corr <- function(y, y_pred) {
  sum(y*y_pred) / (sqrt(sum(y**2)) * sqrt(sum(y_pred**2)))
}

### *** manual_grouping
manual_grouping <- function(x,grpname,expr1){
  grp <- read.table(x, header=TRUE, stringsAsFactors = FALSE)
  rownames(grp) <- grp$vcfid
  return(grp[names(expr1),grpname])
}

### *** pfac_calc
pfac_calc <- function(pfac){
  CMCx2=4.69250216
  ## real max prior
  x2=4.75039291
  y2=0.7
  x1=1.9*x2/CMCx2
  y1=0.7
  ## get penalty factors according to rescaling
  for (j in 1:length(pfac)) {
    pfac[j]=Re((1+y2-2*y1)/((x2-2*x1)^2)*(2*x1^2+pfac[j]*
                                          (x2-2*x1)-2*x1*
                                          (x1^2+pfac[j]*
                                           (x2-2*x1))^0.5)+
               (y1-1)/(x2-2*x1)*(-2*x1+2*(x1^2+pfac[j]*(x2-2*x1))^0.5)+1)
  }
  return(pfac)
}

### ** nested cv for evaluation
nested_cv_elastic_net_perf <- function(x,
                                       y,
                                       n_samples,
                                       n_train_test_folds,
                                       n_k_folds,
                                       alpha,
                                       pfac,
                                       grouping.file) {
  # Gets performance estimates for k-fold cross-validated elastic-net models.
  # Splits data into n_train_test_folds disjoint folds, roughly equal in size,
  # and for each fold, calculates a n_k_folds cross-validated elastic net model. Lambda parameter is
  # cross validated. Then get performance measures for how the model predicts on the hold-out
  # fold. Get the coefficient of determination, R^2, and a p-value, where the null hypothesis
  # is there is no correlation between prediction and observed.
  #
  # The mean and standard deviation of R^2 over all folds is then reported, and the p-values
                                        # are combined using Fisher's method.
  rmse_folds <- rep(0,n_train_test_folds) #@
  R2_folds <- rep(0, n_train_test_folds)
  corr_folds <- rep(0, n_train_test_folds)
  zscore_folds <- rep(0, n_train_test_folds)
  pval_folds <- rep(0, n_train_test_folds)
  ## Outer-loop split into training and test set.
  if(is.na(grouping.file)){
    train_test_fold_ids <- generate_fold_ids(n_samples, n_folds=n_train_test_folds)
  } else {
    train_test_fold_ids <- manual_grouping(grouping.file,"fivefold",y) #@
  }
  for (test_fold in 1:n_train_test_folds) {
    train_idxs <- which(train_test_fold_ids != test_fold)
    test_idxs <- which(train_test_fold_ids == test_fold)
    x_train <- x[train_idxs, ]
    y_train <- y[train_idxs]
    x_test <- x[test_idxs, ]
    y_test <- y[test_idxs]
    ## Inner-loop - split up training set for cross-validation to choose lambda.
    cv_fold_ids <- generate_fold_ids(length(y_train), n_k_folds)
    y_pred <- tryCatch({
      # Fit model with training data.
      fit <- cv.glmnet(x_train,
                       y_train,
                       penalty.factor=pfac,
                       nfolds = n_k_folds,
                       alpha = alpha,
                       type.measure='mse',
                       foldid = cv_fold_ids)
      # Predict test data using model that had minimal mean-squared error in cross validation.
      predict(fit, x_test, s = 'lambda.min')},
      # if the elastic-net model did not converge, predict the mean of the y_train (same as all non-intercept coef=0)
      error = function(cond) rep(mean(y_train), length(y_test)))
    R2_folds[test_fold] <- calc_R2(y_test, y_pred)
    # Get p-value for correlation test between predicted y and actual y.
    # If there was no model, y_pred will have var=0, so cor.test will yield NA.
    # In that case, give a random number from uniform distribution, which is what would
    ## usually happen under the null.
    ##wen's rmse part
    res <- summary(lm(y_test~y_pred)) #@
    rmse_folds[test_fold] <- sqrt(1-res$r.squared)*(res$sigma) #@
    corr_folds[test_fold] <- ifelse(sd(y_pred) != 0, cor(y_pred, y_test), 0)
    zscore_folds[test_fold] <- atanh(corr_folds[test_fold])*sqrt(length(y_test) - 3) # Fisher transformation 
    pval_folds[test_fold] <- ifelse(sd(y_pred) != 0, cor.test(y_pred, y_test)$p.value, runif(1))
  }
  rmse_avg <- mean(rmse_folds) #@
  R2_avg <- mean(R2_folds)
  R2_sd <- sd(R2_folds)
  rho_avg <- mean(corr_folds)
  rho_se <- sd(corr_folds)
  rho_avg_squared <- rho_avg**2
  # Stouffer's method for combining z scores.
  zscore_est <- sum(zscore_folds) / sqrt(n_train_test_folds)
  zscore_pval <- 2*pnorm(abs(zscore_est), lower.tail = FALSE)
  # Fisher's method for combining p-values: https://en.wikipedia.org/wiki/Fisher%27s_method
  pval_est <- pchisq(-2 * sum(log(pval_folds)), 2*n_train_test_folds, lower.tail = F)
  list(rmse_avg = rmse_avg, #@
       R2_avg=R2_avg,
       R2_sd=R2_sd,
       pval_est=pval_est,
       rho_avg=rho_avg,
       rho_se=rho_se,
       rho_zscore=zscore_est,
       rho_avg_squared=rho_avg_squared,
       zscore_pval=zscore_pval)
}

### * others
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
priordf <- fread(priordf.file) ##columns rsid and prior must

genesmodel1 <- list()
genesmodel2 <- list()
for(i in 1:nrow(genes)){
  gene <- genes[i]
  cat("training gene:",gene$gene,"\n")
  gene.ranges <- with(gene, GRanges(seqnames = chr, IRanges(start=start,end=end), gene=gene))
  snpannot.sub <- subsetByOverlaps(snpannot,gene.ranges)
  if(nrow(as.data.frame(snpannot.sub)) < 2){
    r <- data.table(gene=gene$gene,error = "genos is null")
    gene.model <- list(r,r)
  } else {
    index <- as.data.frame(snpannot.sub)$index
    tempgds <- paste0(".temp.",gene$gene,basename(output1),".gds")
    gdsSubset(gds.file,tempgds, snp.include = index)
    mgds <- snpgdsOpen(tempgds)
    genos <- read.gdsn(index.gdsn(mgds, "genotype"))
    ##snpids <- read.gdsn(index.gdsn(mgds,"snp.rs.id"))
    snpids <- as.data.frame(snpannot.sub)$rsid
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
    expr1 <- as.vector(scale(expr1,center=TRUE,scale=TRUE))
    names(expr1) <- expr$sample
    genos <- genos[rownames(genos)%in% names(expr1),]
    expr1 <- expr1[rownames(genos)]
    ## new additions<<<
    ##check functions manual_grouping, pfac_calc
    ##check modifications in nested elastic cv function (lines with #@)
    ##check last argument for prior file should have columns rsid, prior 
    priordf.sub <- data.table(rsid=rsids)
    priordf.sub <- merge(priordf.sub,priordf, by="rsid",all.x=TRUE,
                         sort=FALSE)
    priordf.sub$prior <- ifelse(is.na(priordf.sub$prior),0,priordf.sub$prior)
    priors <- priordf.sub$prior
    names(priors) <- priordf.sub$rsid
    priors <- priors[rsids] ##sort in the same order as genos column
    pfac <- pfac_calc(priors)
    ### new additions >>>
    gene.model <- tryCatch(cvElastic(gene=expr1,
                                     snp=genos,geneid=gene$gene,
                                     pfac=pfac,
                                     grouping.file = grouping.file),
                           error=function(e) {
                             r <- data.table(gene=gene$gene,error=paste0(e))
                             return(list(r,r))
                           } )
    if(!all(is.na(gene.model[[1]]$rsid))){
      gene.model[[1]] <- merge(gene.model[[1]],alleles,by="rsid")
      gene.model[[1]] <- merge(gene.model[[1]],priordf.sub,by="rsid",
                               all.x=TRUE,sort=FALSE)
    } 
  }
  genesmodel1[[i]] <- gene.model[[1]]
  genesmodel2[[i]] <- gene.model[[2]]
}

genesmodel1 <- do.call(bind_rows,genesmodel1)
genesmodel2 <- do.call(bind_rows,genesmodel2)

fwrite(genesmodel1,output1,na="NA")
fwrite(genesmodel2,output2,na="NA")

### * org mode specific
### Local Variables:
### eval: (orgstruct-mode 1)
### orgstruct-heading-prefix-regexp: "### "
### End:
