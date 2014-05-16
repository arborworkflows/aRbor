library(phytools)
library(diversitree)
library(dplyr)
setwd("~/Documents/arbor/aRbor/R")
source("aceArbor.R")
source("makeMkConstraints.R")
source("treeplyr.R")

tree<-pbtree(n=100, scale=1)
Q<-matrix(c(-1,1,1,-1),2,2)
rownames(Q)<-colnames(Q)<-1:2
x<-sim.history(tree,Q)$states
y<-setNames(as.numeric(x),names(x))

ydf <- as.data.frame(y)
ytd<-make.treedata(tree, ydf)

aceArbor(td=ytd, charType="discrete")

aceArbor(td=ytd, colID=2, charType="discrete", discreteModelType="SYM")
aceArbor(ytd, charType="discrete", discreteModelType="ARD") # there is an error here

#aceArbor(tree, as.factor(y), charType="discrete")

#aceArbor(tree, y-1, charType="discrete")


pets<-c("cat", "dog")[y]
names(pets)<-names(y)
petdf<-as.data.frame(pets)
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