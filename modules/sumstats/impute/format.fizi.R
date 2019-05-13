
args <- commandArgs(trailingOnly=TRUE)

if(is.na(args[4])){
  stop("require 4 arguments")
}
sumstats.file <- args[1]
frq.file <- args[2]
out.file <- args[3]
N <- as.numeric(args[4])

library(data.table)

sumstats <- fread(sumstats.file)
sumstats <- sumstats[!P>1] ##Temporary fix
##info filter
sumstats <- sumstats[R2.BLUP>=0.6] ##later can filter for 0.8, if required.

frq <- fread(frq.file)
sumstats <- merge(sumstats,frq[,c("SNP","MAF")],by="SNP")
names(sumstats) <- c("SNP","CHR","BP","A1","A2","TYPE","Z","INFO","P","MAF")
##add N
sumstats$N <- N

##add SE
sumstats$SE <- with(sumstats, 1/sqrt(2 * MAF * (1-MAF) *N))
sumstats$BETA <- sumstats$Z * sumstats$SE
fwrite(sumstats,out.file, sep="\t",na="NA",quote=FALSE)


