
args <- commandArgs(trailingOnly = TRUE)

geneout.file <- args[1]
gwas.file <- args[2]
predixcan.file <- paste0(gwas.file,"_PrediXcan")
epixcan.file <- paste0(gwas.file,"_EpiXcan")
out.file <- args[3]

library(data.table)

dfm <- fread(geneout.file)
dfm$bonp <- p.adjust(dfm$P,method="bonferroni")
dfm$bhp <- p.adjust(dfm$P,method="BH")

predixgenes <- readLines(predixcan.file)
epixgenes <- readLines(epixcan.file)

dfm.bhp <- dfm[bhp<0.05]
dfm.bonp <- dfm[bonp<0.05]

result <- data.table(
  trait=basename(gwas.file),
  predixcan_inMagma_bhp=table(predixgenes %in% dfm.bhp$GENE)["TRUE"],
  predixcan_NotinMagma_bhp=table(predixgenes %in% dfm.bhp$GENE)["FALSE"],
  predixcan_inMagma_bonp=table(predixgenes %in% dfm.bonp$GENE)["TRUE"],
  predixcan_NotinMagma_bonp=table(predixgenes %in% dfm.bonp$GENE)["FALSE"],
  epixcan_inMagma_bhp = table(epixgenes %in% dfm.bhp$GENE)["TRUE"],
  epixcan_NotinMagma_bhp = table(epixgenes %in% dfm.bhp$GENE)["FALSE"],
  epixcan_inMagma_bonp = table(epixgenes %in% dfm.bonp$GENE)["TRUE"],
  epixcan_NotinMagma_bonp = table(epixgenes %in% dfm.bonp$GENE)["FALSE"]
)

fwrite(result,out.file,sep="\t",na="NA",quote=FALSE)
