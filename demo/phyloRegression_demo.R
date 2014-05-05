library(ape)

tree<-pbtree(n=100, scale=1)
mm<-cbind(c(1, -0.5), c(-0.5, 1))
z<-sim.char(tree, mm, model="BM")[,,1]

x<-z[,1]
y<-z[,2]
names(x)<-rownames(z)
names(y)<-rownames(z)
phyloRegression(tree, x, y)