
library(data.table)
library(ggplot2)
library(reshape)

m <- fread("merged.out")
m1 <- reshape(m,idvar="trait",timevar="method",direction="wide")

ggplot(m1,aes(novel.PrediXcan,novel.EpiXcan))+
  geom_point() + geom_abline(intercept = 0, linetype = "dotted") +
  xlab("PrediXcan") + ylab("EpiXcan") + ggtitle("Novel genes identified")
