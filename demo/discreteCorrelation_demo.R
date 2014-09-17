#library(phytools)
#library(diversitree)


tree<-pbtree(n=100, scale=1)
Q<-matrix(c(-1,1,1,-1),2,2)
rownames(Q)<-colnames(Q)<-1:2
x<-sim.history(tree,Q)$states
charA<-setNames(as.numeric(x),names(x))

x<-sim.history(tree,Q)$states
charB<-setNames(as.numeric(x),names(x))

discreteCorrelation(tree, charA, charB)
discreteCorrelation(tree, charA, charB, modelType="ARD")