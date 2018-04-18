#!/bin/env Rscript

### * arguments
args <- commandArgs(trailingOnly = TRUE)
expr.file <- args[1]
covar.file <- args[2]
peerfactors.file <- args[3]
bedout <- args[4]
covarout <- args[5]

### * library
library(data.table)
library("biomaRt")

### * read files
expr <- readRDS(expr.file)
covar <- readRDS(covar.file)
peerfactors <- readRDS(peerfactors.file)

### * get co ordintates and create bed file
genes <- colnames(expr)
genes <- gsub("\\..*$","",genes)
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
##ensembl <- useMart("grch37",dataset="hsapiens_gene_ensembl")

my_refseq <- getBM(attributes=c("ensembl_gene_id","start_position","end_position","chromosome_name","strand"),
                   filters = c("ensembl_gene_id"),
                   values =genes,
                   mart = grch37)
c22 <- my_refseq[my_refseq$chromosome_name == "22",]
c22 <- data.table(c22)
c22$strand <- ifelse(c22$strand == "-1", "-","+")
names(c22) <- c("pid","start","end","Chr","strand")

##change the names of expr
colnames(expr) <- genes
##subset the expression
expr <- expr[,c22$pid]
##transpose
expr <- t(expr)
expr <- as.data.frame(expr)
expr$pid <- rownames(expr)
##merge with c22
expr1 <- merge(c22,expr,by="pid")
##rearrange the columns
expr1 <- data.table(expr1)
expr1$gid <- expr1$pid
f6cols <- c("Chr","start","end","pid","gid","strand")
genenames <- setdiff(names(expr),f6cols)
expr1 <- expr1[,c(f6cols,genenames),with=FALSE]
names(expr1)[1] <- "#Chr"
fwrite(expr1,bedout, sep="\t",na="NA")

### * prepare peer added covar file

## format peer file
peerfactors <- as.data.frame(peerfactors)
peerfactors$id <- rownames(peerfactors)
##rearrange columns
peerfactors <- peerfactors[,c(11,1:10)]

##format covar file
covar <- as.data.frame(covar)
covar$id <- rownames(covar)
##rearrange columns
covar <- covar[,c(15,1:14)]
## merge peer and covar
df <- merge(covar,peerfactors,by="id")
##create a list
peercovars.list <- list()
peercovars.list[[1]] <- covar
for(i in 2:11) peercovars.list[[i]] <- df[,c(names(covar),names(peerfactors)[2:(i)])]

##transpose and write the files in the list
for(i in 1:length(peercovars.list)){
  a <- peercovars.list[[i]]
  a <- t(a)
  a <- as.data.frame(a)
  fwrite(a,paste0(covarout,".",i,".pcovar"),row.names = T,col.names = F,sep="\t",na="NA")
}




### * org mode sepcific
### Local Variables:
### eval: (orgstruct-mode 1)
### orgstruct-heading-prefix-regexp: "### "
### End:

