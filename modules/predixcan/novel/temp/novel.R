#!/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

##argument
clump.file <- args[1]
##e.genes.file <- paste0(dirname(clump.file),"/EpiXcan_",basename(clump.file))
e.genes.file <- paste0(clump.file,"_EpiXcan")
##p.genes.file <- paste0(dirname(clump.file),"/PrediXcan_",basename(clump.file))
p.genes.file <- paste0(clump.file,"_PrediXcan")
annot.file <- args[2]
out.file <- args[3]
window <- 1000000

if(FALSE){
  clump.file <- "PGC2_SCZ.gwas.sumstats.clumped.formatted"
  e.genes.file <- "Sig_genes_for_traits/EpiXcan_PGC2_SCZ"
  p.genes.file <- "Sig_genes_for_traits/PrediXcan_PGC2_SCZ"
  annot.file <- "gencode.v27.build37.txt"
}

getcounts <- function(genes.file,gname){
  library(data.table)
  library(GenomicRanges)
  genes <- readLines(genes.file)
  gwas <- readLines(clump.file, n=1)
  if(gwas =="No Clumps"){
    message("No Clumps in the GWAS file")
    gwas <- data.table(SNP="rs0000",
                       BP=0,
                       CHR=1)
  } else {
    gwas <- fread(clump.file)[,c("SNP","BP","CHR")]
  }
  nlocus <- nrow(gwas) ##
  gcode <- fread(annot.file)

  table(genes %in% gcode$gene_id)
  genes.gcode <- gcode[gene_id %in% genes]

  genes.r <- with(genes.gcode,
                  GRanges(seqnames=seqnames,
                          IRanges(start=start,
                                  end=end),
                          strand=strand,
                          geneID = gene_id,
                          gene = gene_name
                          ))

  ##covert gwas to ranges
  gwas.r <- with(gwas, GRanges(seqnames=paste0("chr",CHR),
                               IRanges(start=BP-window,
                                       end=BP+window)))

  genes.r1 <- subsetByOverlaps(genes.r,gwas.r)
  genes.r1 <- as.data.frame(genes.r1)
  known <- nrow(genes.r1)
  novel <- nrow(genes.gcode) - known
  known.genes  <- paste(genes.r1$geneID,sep=":")

  return(data.table(
    trait = gsub(".clumped.formatted","",clump.file),
    method = gname,
    ngenes = nrow(genes.gcode),
    known = known,
    novel = novel,
    known.genes = known.genes
  ))
}

e.res <- getcounts(genes.file=e.genes.file,gname="EpiXcan")
p.res <- getcounts(genes.file=p.genes.file,gname="PrediXcan")

res.df <- rbind(e.res,p.res)
fwrite(res.df,out.file)
