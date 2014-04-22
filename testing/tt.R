library(expm)
qq<-cbind(c(-1, 1), c(1, -1))

library(diversitree)

exTree<-read.tree(text="((A:1.0,B:1.0):4.0,(C:2.0,(D:1.5,E:1.5):0.5):3.0);")
exData<-c(1, 2, 2, 1, 1)
names(exData)<-c("A", "B","C", "D", "E")

lik<-make.mkn(exTree, exData, k=2)
lik<-constrain(lik, q12~q21)


lik(1, root=ROOT.FLAT) # should be -3.482415
lik(0.1, root=ROOT.FLAT) # should be -5.006
		
fit<-find.mle(lik, setNames(1, argnames(lik)), root=ROOT.FLAT)

# make sure root is being passed through
# just so that I know what I'm doing

library(geiger)
library(TreeSim)

tt<-sim.bd.taxa.age(100, 1, 1, 0, age=10)[[1]]
plot(tt)

cc<-sim.character(tt, c(.1, .1), x0=0, model="mk2")+1

lik<-make.mkn(tt, cc, k=2)
lik<-constrain(lik, q12~q21)

lik(0.1)
lik(0.1, root=ROOT.FLAT)

fit1<-find.mle(lik, setNames(1, argnames(lik)))
fit2<-find.mle(lik, setNames(1, argnames(lik)), root=ROOT.FLAT)

# these are different
fit1$lnL
lik(fit1$par)
lik(fit1$par, root=ROOT.FLAT)

fit2$lnL
lik(fit2$par)
lik(fit2$par, root=ROOT.FLAT)

# these results are as expected

