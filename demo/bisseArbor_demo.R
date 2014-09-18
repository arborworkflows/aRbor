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

aa<-bisseArbor(td=ytd)
