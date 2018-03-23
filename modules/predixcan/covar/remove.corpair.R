##replace covar names with dummy names
newcovar <- covar[,-c("id"),with=FALSE]
covar.names <- names(newcovar)
covar.names <- data.frame(covar.names)
covar.names$dummy <- paste0("cov",1:nrow(covar.names))
names(newcovar) <- covar.names$dummy

##pairwise correlation
fmula <- as.formula(paste0("~",paste0(names(newcovar),collapse = "+")))

c <- tryCatch(canCorPairs(fmula,newcovar),
              error=function(e){
                r <- gsub("Error.*.analyzed:","",paste0(e))
                r <- gsub("\n","",r)
                r <- gsub(" ","",r)
                r <- unlist(strsplit(r,split=","))
                message("removing ",r)
                newcovar <- newcovar[,-r,with=F]
                fmula <- as.formula(paste0("~",paste0(names(newcovar),collapse = "+")))
                return(canCorPairs(fmula,newcovar))
              })

##remap the names
newcovar.names <- data.frame(dummy=rownames(c))
newcovar1.names <- merge(newcovar.names,covar.names,by="dummy",sort=FALSE,all.x=TRUE)
if(!all(newcovar1.names$dummy==newcovar.names$dummy)) stop("error in merge")
c1 <- c
rownames(c1) <- newcovar1.names$covar.names
colnames(c1) <- newcovar1.names$covar.names

pdf("pairwise.pdf")
plotCorrMatrix(c1)
dev.off()

##melt the c
cmelt <- melt(c)
cmelt <- data.table(cmelt)
cmelt$X1 <- as.character(cmelt$X1)
cmelt$X2 <- as.character(cmelt$X2)
cmelt <- cmelt[!X1==X2]

##remove one of each correlated pairs
toremove <- c()
tokeep <- c("Dx","Institution","gPC1","gPC2","gPC3","gPC4")

for(i in 1:nrow(cmelt)){
  cpairs <- c(cmelt$X1[i],cmelt$X2[i])
  if(any(cpairs %in% toremove)) {
    pass
  }
}
