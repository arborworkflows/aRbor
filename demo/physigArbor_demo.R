#library(phytools)
#library(diversitree)
#library(geiger)
#library(picante)

tree<-pbtree(n=100, scale=1)
Q<-matrix(c(-1,1,1,-1),2,2)
rownames(Q)<-colnames(Q)<-1:2
x<-sim.history(tree,Q)$states
y<-setNames(as.numeric(x),names(x))

td<-make.treedata(tree, y)

physigArbor(td, charType="discrete")
physigArbor(td, charType="discrete", signalTest="garbageTest")
#physigArbor(tree, y, charType="discrete", signalTest="Blomberg") # would return error

z<-sim.char(tree, 1, model="BM")[,,1]
td2<-make.treedata(tree, z)

physigArbor(td2, charType="continuous")
physigArbor(td2, charType="continuous", signalTest="garbageTest")
physigArbor(td2, charType="continuous", signalTest="Blomberg")
