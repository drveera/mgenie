library(data.table)

args <- commandArgs(trailingOnly = TRUE)

gwas.file <- args[1]
out.file <- args[2]

gwas <- fread(gwas.file)
gwas <- gwas[CHR!=22]
gwas$fdr <- p.adjust(gwas$P, method="BH")
gwas.sig <- gwas[P<0.05]
gwas.sig1 <- gwas[P<0.01]
gwas.fdr <- gwas[fdr<0.05]

res <- data.table(
  pheno=out.file,
  fdr=nrow(gwas.fdr),
  nominal=nrow(gwas.sig),
  nominal1=nrow(gwas.sig1)
)

fwrite(res,out.file,sep=",",na="NA",col.names = FALSE)
