#!/bin/env Rscript

### * arguments
args <- commandArgs(trailingOnly = TRUE)
residuals.file <- args[1]
bed.file <- args[2]
bedout <- args[3]


### * library
library(data.table)
library("biomaRt")

### * read files
expr <- readRDS(residuals.file)
##genes in rows and ids in columns, so transpose it
expr <- t(expr)

### * read the input bed file
c22 <- fread(bed.file)

##change the names of expr
genes <- colnames(expr)
genes <- gsub("\\..*$","",genes)
colnames(expr) <- genes

##subset the expression
genescommon <- intersect(c22$pid,genes)
expr <- expr[,genescommon]
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


### * org mode sepcific
### Local Variables:
### eval: (orgstruct-mode 1)
### orgstruct-heading-prefix-regexp: "### "
### End:

