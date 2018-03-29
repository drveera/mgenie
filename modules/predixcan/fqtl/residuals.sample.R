#!/bin/env Rscript

### * libraries
library(limma)
library(data.table)
library(edgeR)

### * input files
args <- commandArgs(trailingOnly = TRUE)
expr.file <- args[1]
covar.file <- args[2]
out.file <- args[3]


### * functions
### ** create model matrix
createDm <- function(x){
  x <- x[,-c("id"),with=FALSE]
  x.classes <- summary.default(x)[,"Mode"]
  x.classes.numeric <- x.classes[x.classes=="numeric"]
  x.classes.char <- x.classes[x.classes=="character"]
  x1 <- x[,names(x.classes.char),with=FALSE]
  x2 <- x[,names(x.classes.numeric),with=FALSE]
  fm <- as.formula(paste0("~",paste0(names(x1),collapse = "+")))
  x1.dm <- model.matrix(fm,data=x1)
  print(dim(x1.dm))
  print(dim(x))
  return(cbind(x1.dm,as.matrix(x2)))
}

### ** writeandread
writeandread <- function(x){
  fwrite(x,".temp")
  return(fread(".temp"))
}

### ** calcResiduals
calcResiduals <- function(geneBySampleValues, samplesByCovariates, varsToAddBackIn=NULL, sampleWeights=NULL) {
  #################################################################################
  # Use the calcResiduals() code after this section in a for loop:
  #################################################################################
  if (is.matrix(sampleWeights)) {
    residualizedMat = matrix(NA, nrow=nrow(geneBySampleValues), ncol=ncol(geneBySampleValues), dimnames=dimnames(geneBySampleValues))
    for (gInd in 1:nrow(geneBySampleValues)) {
      gRow = calcResiduals(geneBySampleValues[gInd, , drop=FALSE], samplesByCovariates, varsToAddBackIn, sampleWeights[gInd, ])
      residualizedMat[gInd, ] = gRow
    }
    return(residualizedMat)
  }
  #################################################################################
  
  #result.lm = lsfit(x=samplesByCovariates, y=t(geneBySampleValues), wt=sampleWeights, intercept=FALSE)
  
  # Formula of "y ~ 0 + x" means no intercept:
  result.lm = lm(t(geneBySampleValues) ~ 0 + samplesByCovariates, weights=sampleWeights)
  covarNames = colnames(samplesByCovariates)
  
  coef = result.lm$coefficients
  isMatrixForm = is.matrix(coef)
  if (isMatrixForm) {
    rownames(coef) = covarNames
  }
  else {
    names(coef) = covarNames
  }
  
  allVarsToAddBack = '(Intercept)'
  if (!is.null(varsToAddBackIn)) {
    allVarsToAddBack = c(allVarsToAddBack, varsToAddBackIn)
  }
  allVarsToAddBack = intersect(allVarsToAddBack, covarNames)
  
  residualizedMat = result.lm$residuals
  for (v in allVarsToAddBack) {
    if (isMatrixForm) {
      multCoef = coef[v, , drop=FALSE]
    }
    else {
      multCoef = coef[v]
    }
    residualizedMat = residualizedMat + samplesByCovariates[, v, drop=FALSE] %*% multCoef
  }
  
  residualizedMat = t(residualizedMat)
  rownames(residualizedMat) = rownames(geneBySampleValues)
  colnames(residualizedMat) = colnames(geneBySampleValues)
  
  return(residualizedMat)
}

### * read files
### ** read expression matrix
dge <- readRDS(expr.file)

### ** read and process covar
covar <- fread(covar.file)
## convert to matrix
covar.ids <- covar$id
covar <- covar[,-c("id")]
covar <- as.matrix(covar)
rownames(covar) <- covar.ids
## make sure samples in covar is in same order as expr file
print(head(covar))
covar <- covar[colnames(dge),]

### * analysis
### ** run voom
expr.voom <- voomWithQualityWeights(dge,covar,plot=TRUE,normalize.method = "none")

### ** calculate residuals
expr.residuals <- calcResiduals(expr.voom$E,covar, sampleWeights = expr.voom$weights)

### output
saveRDS(expr.residuals, out.file)

### * org mode specific
### Local Variables:
### eval: (orgstruct-mode 1)
### orgstruct-heading-prefix-regexp: "### "
### End:
