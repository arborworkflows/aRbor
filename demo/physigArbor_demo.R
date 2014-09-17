#library(phytools)
#library(diversitree)
#library(geiger)
#library(picante)

tree<-pbtree(n=100, scale=1)
Q<-matrix(c(-1,1,1,-1),2,2)
rownames(Q)<-colnames(Q)<-1:2
x<-sim.history(tree,Q)$states
y<-setNames(as.numeric(x),names(x))

physigArbor(tree, y, charType="discrete")
physigArbor(tree, y, charType="discrete", signalTest="garbageTest")
physigArbor(tree, y, charType="discrete", signalTest="Blomberg")

z<-sim.char(tree, 1, model="BM")[,,1]
physigArbor(tree, z, charType="continuous")
physigArbor(tree, z, charType="continuous", signalTest="garbageTest")
physigArbor(tree, z, charType="continuous", signalTest="Blomberg")
