library(phytools)

tree<-pbtree(n=100, scale=1)
Q<-matrix(c(-2,1,1,1,-2,1,1,1,-2),3,3)
rownames(Q)<-colnames(Q)<-1:3
x<-sim.history(tree,Q)$states

# ASR using rerootingMethod
XX<-rerootingMethod(tree,x,model="ER")
# ASRs in XX$marginal.anc

# ASR using diversitree
y<-setNames(as.numeric(x),names(x))
# (we needed to convert to numeric)
lik<-make.mkn(tree,y,k=3)
# constrain to ER model
lik<-constrain(lik,q13~q12,q21~q12,q23~q12,q31~q12,q32~q12)
fit<-find.mle(lik,setNames(1,argnames(lik)))
ZZ<-t(asr.marginal(lik,coef(fit)))
dimnames(ZZ)<-dimnames(XX$marginal.anc)

ZZ<-t(asr.joint(lik,coef(fit)))

fd<-fitDiscrete(tree, y)

# Compare
XX$loglik
fit$lnLik
fd$opt$lnL

XX$Q
fit$par
fd$opt$q12

fd<-fitDiscrete(tree, y)

tree<-pbtree(n=100, scale=1)
Q<-matrix(c(-1,1,1,-1),2,2)
rownames(Q)<-colnames(Q)<-1:2
x<-sim.history(tree,Q)$states
y<-setNames(as.numeric(x),names(x))

aceArbor(tree, y, charType="discrete")