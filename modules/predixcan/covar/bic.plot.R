#!/bin/env Rscript

### * arguments
args <- commandArgs(trailingOnly = TRUE)
bic.file <- args[1]
out.file1 <- args[2]
out.file2 <- args[3]

### * libraries
library(data.table)
library(ggplot2)

### * read files
bic.df <- fread(bic.file)

### * process data
model.df <- data.frame(table(bic.df$bestModel))
model.df <- data.table(model.df)
model.df0 <- model.df[order(Freq,decreasing = T)]
model.df <- model.df0[1:20]
model.df$Var1 <- factor(model.df$Var1, levels = model.df$Var1)

### * individuals covars
models <- as.character(model.df0$Var1)
models <- gsub(" ","",models)
temp <- unlist(strsplit(models, split="\\+"))
temp.df <- data.table(data.frame(table(temp)))
temp.df <- temp.df[order(Freq, decreasing = TRUE)]
temp.df$temp <- factor(temp.df$temp, levels=temp.df$temp)

### * individual covar plot
p2 <- ggplot(temp.df, aes(temp,Freq)) + geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p1 <- ggplot(model.df, aes(Var1, Freq)) + geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot=p1,filename=out.file1)
ggsave(plot=p2,filename=out.file2)
print("done")


### * org mode specific
### Local Variables:
### eval: (orgstruct-mode 1)
### orgstruct-heading-prefix-regexp: "### "
### End:
