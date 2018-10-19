
args <- commandArgs(trailingOnly = TRUE)


bed.file <- args[1]
out.file <- args[2]


library(data.table)
##format bed file
ebed <- fread(paste0("zcat ",bed.file))
names(ebed)[1] <- "chr"
ebed <- ebed[chr==22]
ebed <- ebed[,7:ncol(ebed)]
ebed <- t(ebed)
famids <- rownames(ebed)
ebed <- data.table(ebed)
pheno <- cbind(famids,famids,ebed)

for(i in 1:ncol(ebed)){
  tpheno <- paste0(out.file,".",i)
  fwrite(pheno[,c(1,2,(2+i)),with=F],tpheno,
         sep=" ",na="NA",col.names = FALSE)
}





