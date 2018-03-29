#!/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

count.file <- args[1]
outfile <- args[2]

df <- read.table(count.file, stringsAsFactors=FALSE)
names(df) <- c("count","filename")
df$filename <- basename(df$filename)
df$filename <- gsub("^.*RDS.","",df$filename)
df$filename <- gsub(".nominalout","",df$filename)
df$filename <- as.numeric(df$filename)
df <- df[order(df$filename),]
df$filename <- paste0(df$filename,"_peer_factors")
df$filename <- factor(df$filename, levels=df$filename)

library(ggplot2)

p1 <- ggplot(df, aes(filename,count, group=1)) + geom_point() +
  geom_line() + theme(axis.text.x=element_text(angle=90)) + xlab("Number of peer factors") + ylab("eQTLs with P < 0.01")

ggsave(plot=p1, filename=outfile)
