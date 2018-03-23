#!/bin/env Rscript

### * arguments 
args <- commandArgs(trailingOnly = TRUE)
expr.file <- args[1]
covar.file <- args[2]
out.file.rds <- args[3]
out.file.covar <- args[4]
genelist.file <- args[5]
thold <- as.numeric(args[6])
out.file.voom <- args[7]
cutoff <- as.numeric(args[8])

### * test arguments commented
##expr.file <- "sample.raw.RDS"
##covar.file <- "covar.newPCs.txt"
##out.file.rds <- "testout.rds"
##out.file.covar <- "testout.covar"
##genelist.file <- "testout.genelist"
##thold <- 1
##out.file.voom <- "testout.voom"

### * libraries
library(data.table,warn.conflicts = FALSE)
library(reshape, warn.conflicts = FALSE)
library(limma, warn.conflicts = FALSE, verbose = FALSE)
library(edgeR)
library(ggplot2)
library(plotly, warn.conflicts = FALSE)
library(variancePartition)

### * functions
### ** removepairs: remove one from each pairs
removepairs <- function(c,tokeep,thold){
  cmelt <- melt(c)
  cmelt <- data.table(cmelt)
  cmelt$X1 <- as.character(cmelt$X1)
  cmelt$X2 <- as.character(cmelt$X2)
  cmelt <- cmelt[!X1==X2]
  cmelt <- cmelt[value>thold]

  ##remove one of each correlated pairs
  tokeep <- c()
  toremove <- c()
  for(i in 1:nrow(cmelt)){
    cpairs <- c(cmelt$X1[i],cmelt$X2[i])
    cpairs <- cpairs[!cpairs %in% tokeep]
    if(length(cpairs)==0 | any(cpairs %in% toremove)) next
    preserve <- c("RNAFlowcell_Batch","RNAIntronic_Rate","RNALibrary_Batch")
    if (cpairs[1] %in% preserve){
      toremove <- c(toremove, cpairs[2])
      tokeep <- c(tokeep,cpairs[1])
    } else {
      toremove <- c(toremove,cpairs[1])
      tokeep <- c(tokeep, cpairs[2])
    }
  }
  writeLines(tokeep,paste0(out.file.covar,".CorVars.kept"))
  return(toremove)
}

### ** paircor: pairwise correlation
paircor <- function(covar,plname){ 
  newcovar <- covar[,-c("id"),with=FALSE]
  covar.names <- names(newcovar)
  covar.names <- data.frame(covar.names)
  covar.names$dummy <- paste0("cov",1:nrow(covar.names)) ##replace covar names with dummy names
  names(newcovar) <- covar.names$dummy

  fmula <- as.formula(paste0("~",paste0(names(newcovar),collapse = "+")))
  c <- tryCatch(canCorPairs(fmula,newcovar),
                error=function(e){
                  print(e)
                  r <- gsub("Error.*.analyzed:","",paste0(e))
                  r <- gsub("\n","",r)
                  r <- gsub(" ","",r)
                  r <- unlist(strsplit(r,split=","))
                  message("removing ",r)
                  newcovar <- newcovar[,-r,with=F]
                  fmula <- as.formula(paste0("~",paste0(names(newcovar),collapse = "+")))
                  return(canCorPairs(fmula,newcovar))
                })
  ##remap the names
  newcovar.names <- data.frame(dummy=rownames(c))
  newcovar1.names <- merge(newcovar.names,covar.names,by="dummy",sort=FALSE,all.x=TRUE)
  if(!all(newcovar1.names$dummy==newcovar.names$dummy)) stop("error in merge")
  c1 <- c
  rownames(c1) <- newcovar1.names$covar.names
  colnames(c1) <- newcovar1.names$covar.names

  ##pdf("pairwise.pdf")
  pdf(paste0(plname,".preRemoval.pdf"))
  plotCorrMatrix(c1)
  dev.off()

  ## remove one from each pair
  ##  tokeep <- c("Dx","Institution","gPC1","gPC2","gPC3","gPC4") #from function args
  tokeep <- c()
  toremove <- removepairs(c1,tokeep,thold=thold) ##custom function
  writeLines(toremove,paste0(out.file.covar,".Cor.Vars.removed"))
  cnames <- rownames(c1)
  tokeep <- unique(c(tokeep,cnames[!cnames %in% toremove]))
  c1.keep <- c1[tokeep,tokeep]
  pdf(paste0(plname,".postRemoval.pdf"))
  plotCorrMatrix(c1.keep)
  dev.off()
  return(covar[,c("id",tokeep),with=FALSE])
}

### * load files
expr <- readRDS(expr.file)
covar <- fread(covar.file)

### * Normalization
### ** remove all NAs
covar <- covar[complete.cases(covar),]

### ** note down the column names
covar.names <- names(covar)
covar.names <- covar.names[!covar.names %in% c("id")]

### ** check if the samples are in rows
ids.merge <- intersect(rownames(expr),covar$id)
if(length(ids.merge)==0) stop("no samples from covar file match with rownames of expression file")
### ** subset expression
expr <- expr[ids.merge,]

### ** transpose the expression file
exprt <- t(expr)


### ** run normalization
### *** dge
dge <- DGEList(exprt)
### ** filter out low counts
genestokeep <- rowSums(cpm(dge) > 1) >= (ncol(dge)*cutoff)
cat("using cutoff",cutoff,"\n")
cat("total number of final genes ",table(genestokeep)["TRUE"])
dge <- dge[genestokeep,]

### *** tmm
dge <- calcNormFactors(dge,method="TMM")

### *** voom and plot
pdf(paste0(out.file.covar,".voom.pdf"))
expr.voomed <- voom(dge,plot=TRUE)
dev.off()
### *** get voomed expression
expr <- expr.voomed$E
### *** transpose the expression
expr <- t(expr) ##this is used for pca and save for next step
### ** export files
saveRDS(expr,out.file.rds)
saveRDS(expr.voomed,out.file.voom)
writeLines(colnames(expr),genelist.file)
### * Expression pca
### ** run pca and get significant pcs
exprt.pca <- prcomp(expr,retx=TRUE)
exprt.pcs <- exprt.pca$x
##p1 <- ggplot(as.data.frame(exprt.pcs),aes(PC1,PC2)) + geom_point()
##ggsave(plot=p1,filename=paste0(out.file,".PCs1and2.pdf"))
pca.importance <- summary(exprt.pca)$importance
pca.importance <- t(pca.importance)
pca.ids <- rownames(pca.importance)
pca.importance <- data.table(pca.importance)
pca.importance$pcs <- pca.ids
names(pca.importance) <- gsub(" ",".",names(pca.importance))
pca.importance.sig <- pca.importance[Proportion.of.Variance>=0.1]
if(nrow(pca.importance.sig)==0){
  pca.importance.sig <- pca.importance[Proportion.of.Variance>=0.05]
  if(nrow(pca.importance.sig)==0){
    stop("no pcs had variance >= 0.05")
  }
}
##subset pcs
goodpcs <- pca.importance.sig$pcs
exprt.pcs <- as.matrix(exprt.pcs[,goodpcs])
colnames(exprt.pcs) <- goodpcs
pcids <- rownames(exprt.pcs)
exprt.pcs <- data.table(exprt.pcs)
exprt.pcs$id <- pcids

### * PCA association with covariates
### ** merge with covars
covar.noPC <- covar
covar <- merge(covar,exprt.pcs,by="id")
### ** assoc
fullres <- list()
for(j in 1:length(goodpcs)){
  pc <- goodpcs[j]
  res <- list()
  for(i in 1:length(covar.names)){
    cov.name <- covar.names[i]
    df <- covar[,c(pc,cov.name),with=FALSE]
    names(df) <- c("pc","cov")
    print(head(df))
    m <- glm("pc~cov",data=df)
    r2 <- fmsb::NagelkerkeR2(m)
    m.coeff <- as.data.frame(summary(m)$coeff)[-1,]
    names(m.coeff)[4] <- "pvalue"
    m.coeff <- m.coeff[m.coeff$pvalue==min(m.coeff$pvalue),][1,]
    pvalue <- m.coeff$pvalue
    svar <- rownames(m.coeff)
    res[[i]] <- data.frame(pc=pc,
                           covar=cov.name,
                           svar=svar,
                           pvalue=pvalue,
                           r2=r2)
  }
  res <- do.call(rbind,res)
  res$pvalue <- ifelse(is.na(res$pvalue),1,res$pvalue)
  fullres[[j]] <- res
}

fullres <- do.call(rbind,fullres)
fullres <- data.table(fullres)

### ** plot the associations commented
##pc.assoc.plot <- ggplot(fullres, aes(r2.R2,-log10(pvalue), colour=pc,text=covar)) + geom_point(size=2)
##pc.assoc.plotly <- ggplotly(pc.assoc.plot)
##htmlwidgets::saveWidget(pc.assoc.plotly,"~/bridge/pc.assoc.html")
##htmlwidgets::saveWidget(pc.assoc.plotly,paste0(getwd(),"/",out.file.covar,".PCA.assoc.html"))
 
### ** select the significant covariates
full.res.sig <- fullres[pvalue<0.05]
full.res.nosig <- fullres[pvalue>0.05]
selected.covariates <- unique(as.character(full.res.sig$covar))
unselected.covariates <- unique(as.character(full.res.nosig$covar))
writeLines(unselected.covariates,paste0(out.file.covar,".removedDueToExprPCS"))

### * Correlation between covariates
sel.covar <- paircor(covar.noPC,out.file.covar)
fwrite(sel.covar,out.file.covar,sep="\t",na="NA")


### * org mode specific
### Local Variables:
### eval: (orgstruct-mode 1)
### orgstruct-heading-prefix-regexp: "### "
### End:
