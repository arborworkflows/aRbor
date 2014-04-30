library(phytools)
library(diversitree)


tree<-pbtree(n=100, scale=1)
Q<-matrix(c(-1,1,1,-1),2,2)
rownames(Q)<-colnames(Q)<-1:2
x<-sim.history(tree,Q)$states
y<-setNames(as.numeric(x),names(x))

aceArbor(tree, y, charType="discrete")
aceArbor(tree, y, charType="discrete", discreteModelType="SYM")
aceArbor(tree, y, charType="discrete", discreteModelType="ARD") # there is an error here

aceArbor(tree, as.factor(y), charType="discrete")

aceArbor(tree, y-1, charType="discrete")


zz<-c("cat", "dog")[y]
names(zz)<-names(y)
aceArbor(tree, zz, charType="discrete")
aceArbor(tree, zz, charType="discrete", aceType="joint")
aceArbor(tree, zz, charType="discrete", aceType="MCMC") # this doesn't look right
 
