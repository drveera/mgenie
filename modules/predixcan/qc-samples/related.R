#!/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

genome.file <- args[1]
fam.file <- args[2]
output.file <- args[3]


library(data.table)

##identify ids to remove
genome <- fread(genome.file,header=TRUE)
genome <- genome[PI_HAT>=0.2]
if(nrow(genome)==0){
  fam <- fread(fam.file)
  fam <- fam[V1=="???"] ##just make the df empty
  fwrite(fam,output.file,sep="\t",na="NA")
} else {
  genome <- genome[IID1 != IID2]
  relpairs <- genome[,c("IID1","IID2"),with=FALSE]
  allids <- c(relpairs$IID1,relpairs$IID2)
  frq <- data.frame(table(allids))
  print(head(frq))
  names(frq) <- c("IID1","IID1.freq")
  frq$IID1 <- as.character(frq$IID1)
  relpairs$IID1 <- as.character(relpairs$IID1)
  relpairs <- merge(relpairs,frq,by="IID1",all.x=TRUE,sort=FALSE)
  names(frq) <- c("IID2","IID2.freq")
  frq$IID2 <- as.character(frq$IID2)
  relpairs$IID2 <- as.character(relpairs$IID2)
  relpairs <- merge(relpairs,frq,by="IID2",all.x=TRUE,sort=FALSE)
  relpairs$IID1 <- as.character(relpairs$IID1)
  relpairs$IID2 <- as.character(relpairs$IID2)

  toremoveids <- function(relpairs){
    toremove <- c()
    for(i in 1:nrow(relpairs)){
      id1 <- relpairs$IID1[i]
      id2 <- relpairs$IID2[i]
      id1.f <- relpairs$IID1.freq[i]
      id2.f <- relpairs$IID2.freq[i]
      id1.present <- id1 %in% toremove
      id2.present <- id2 %in% toremove
      if(!(id1.present|id2.present)){
        if(id1.f > id2.f){
          toremove <- c(toremove,id1)
        } else {
          toremove <- c(toremove,id2)
        }
      }
    }
    return(toremove)
  }

  ##reorder decreasing IID1
  relpairs <- relpairs[order(IID1.freq,decreasing = TRUE)]
  toremove1 <- toremoveids(relpairs)

  ##reorder decreasing IID2
  relpairs <- relpairs[order(IID2.freq,decreasing = TRUE)]
  toremove2 <- toremoveids(relpairs)


  if(length(toremove1)<length(toremove2)){
    toremove <- toremove1
  } else {
    toremove <- toremove2
  }

  ##subset fam
  fam <- fread(fam.file)
  fam <- fam[V2 %in% toremove]

  fwrite(fam,output.file,sep="\t",na="NA")
}
