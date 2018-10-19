args <- commandArgs(trailingOnly = TRUE)

library(data.table)

tfile <- args[1]
out.file <- args[2]

tstat <- fread(tfile)
names(tstat) <- c("pheno","fdr","nominal","nominal2")

pheno <- gsub("^.*RDS.","",tstat$pheno)
peer <- as.numeric(gsub("\\..*$","",pheno))

tstat$peer <- as.numeric(peer)

##fdr plot

peer.split <- split(tstat,by="peer")

res <- lapply(peer.split, function(x){
  fdr <- table(x$fdr>0)["TRUE"]
  fdr <- ifelse(is.na(fdr),0,fdr)
  sig5.mean <- mean(x$nominal,na.rm = TRUE)
  sig5.SD <- sd(x$nominal,na.rm = TRUE)
  sig1.mean <- mean(x$nominal2,na.rm = TRUE)
  sig1.SD <- sd(x$nominal2,na.rm = TRUE)
  fdr.mean <- mean(x$fdr,na.rm=TRUE)
  return(data.table(fdr=fdr,
                    sig5.mean=sig5.mean,
                    sig5.SD=sig5.SD,
                    sig1.mean=sig1.mean,
                    sig1.SD=sig1.SD,
                    fdr.mean=fdr.mean))
})

for(i in 1:length(res)) res[[i]]$peer <- names(res)[i]
res <- do.call(rbind,res)
res$peer <- as.numeric(res$peer)
res$peer <- factor(res$peer,levels=0:max(res$peer))

library(ggplot2)

p1 <- ggplot(res, aes(peer,fdr)) + geom_point() +
  geom_line(aes(group=1)) +
  ggtitle("No. of genes with atleast one trans QTL")
ggsave(plot=p1,filename=paste0(out.file,".fdr.pdf"))


tstat$peer <- factor(tstat$peer, levels=0:max(peer))
p3 <- ggplot(tstat, aes(peer,nominal)) + geom_violin() +
  geom_point(data=res,aes(peer,sig5.mean)) +
  geom_line(data=res,aes(peer,sig5.mean,group=1))+
  ggtitle("No. of genes with P value < 0.05")
ggsave(plot=p3, filename=paste0(out.file,"p0.05.violin.pdf"))

p4 <- ggplot(tstat, aes(peer,nominal2)) + geom_violin() +
  geom_point(data=res,aes(peer,sig1.mean)) +
  geom_line(data=res,aes(peer,sig1.mean,group=1))+
  ggtitle("No. of genes with P value < 0.01")
ggsave(plot=p4,filename=paste0(out.file,"p0.01.violin.pdf"))


p5 <- ggplot(tstat, aes(peer,fdr)) + geom_violin()+
  geom_point(data=res,aes(peer,fdr.mean)) +
  geom_line(data=res,aes(peer,fdr.mean,group=1))+
  ggtitle("No. of genes with FDR P  < 0.01")
ggsave(plot=p5,filename=paste0(out.file,".FDR.violin.pdf"))

