
args  <- commandArgs(trailingOnly = TRUE)

hsq.file  <- args[1]
out.file  <- args[2]

library(data.table)

hsq  <- fread(hsq.file, fill=TRUE)

h2  <- hsq[Source=="V(G)/Vp"]$Variance[1]
se  <-hsq[Source=="V(G)/Vp"]$SE[1]
pvalue  <- hsq[Source=="Pval"]$Variance[1]
n  <- hsq[Source=="n"]$Variance[1]

rdf  <- data.table(filename=basename(hsq.file),
                   h2 = h2,
                   se = se,
                   pvalue = pvalue,
                   n = n)

fwrite(rdf,out.file,na="NA")
