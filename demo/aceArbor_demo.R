#library(phytools)
#library(diversitree)
#library(dplyr)
#setwd("~/Documents/arbor/aRbor/R")
#source("aceArbor.R")
#source("makeMkConstraints.R")
#source("treeplyr.R")

tree<-pbtree(n=100, scale=1)
Q<-matrix(c(-1,1,1,-1),2,2)
rownames(Q)<-colnames(Q)<-1:2
x<-sim.history(tree,Q)$states
y<-setNames(as.numeric(x),names(x))

ydf <- as.data.frame(y)
ytd<-make.treedata(tree, ydf)

aa<-aceArbor(td=ytd, charType="discrete")
plot(aa)

aa<-aceArbor(td=ytd, charType="discrete", aceType="stochastic")
plot(aa[[1]], ytd$phy)

aceArbor(td=ytd, charType="discrete", discreteModelType="SYM")
aceArbor(ytd, charType="discrete", discreteModelType="ARD")



pets<-c("cat", "dog")[y]
names(pets)<-names(y)
x2<-sim.history(tree,Q)$states
y2<-setNames(as.numeric(x),names(x))
pets2<-c("pork", "beef")[y2]

petdf<-as.data.frame(cbind(pets, pets2))
pettd<-make.treedata(tree, petdf)

aceArbor(pettd, charType="discrete")
aceArbor(pettd, charType="discrete", aceType="joint")
aceArbor(pettd, charType="discrete", aceType="MCMC") # this doesn't look right
 


tree<-pbtree(n=50)
x<-fastBM(tree) # simulate using fastBM
td<-make.treedata(tree, x)

x<-as.data.frame(x)
td<-make.treedata(tree, x)

aceArbor(td)