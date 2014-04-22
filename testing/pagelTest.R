pagelTree<-read.tree(text="((s1:1.0,s2:1.0):1.0,s3:2.0);")
pagelData<-c(1, 2, 2)
names(pagelData)<-c("s1", "s2","s3")

lik<-make.mkn(pagelTree, pagelData, k=2)
lik<-constrain(lik, q12~q21)
		
fit<-find.mle(lik, setNames(1, argnames(lik)))

		zz<-asr.marginal(lik, coef(fit))
		
		plot(phy, main="Marginal ASR")
		nodelabels(pie=t(zz), piecol=1:2, cex=.5, frame="circle")
		
		return(t(zz))