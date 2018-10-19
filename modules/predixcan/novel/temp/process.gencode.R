
library(data.table)
library(GenomicRanges)

temp <- rtracklayer::import("gencode.v27lift37.basic.annotation.gtf.gz")
temp <- as.data.frame(temp)
temp <- data.table(temp)
temp <- temp[,c("seqnames","start","end","strand","gene_id","gene_name")]
temp$gene_id <- gsub("\\..*$","",temp$gene_id)
temp <- temp[!duplicated(gene_id)]

fwrite(temp,"gencode.v27.build37.txt")

temp.r <- with(temp,GRanges(seqnames = seqnames,
                            IRanges(start=start,end=end),
                            strand=strand,
                            geneID=gene_id,
                            gene=gene_name))
