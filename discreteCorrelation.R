discreteCorrelation<-function(phy, d1, d2) {
	
	# check for perfect name matching
	# this could be improved a lot
	
	td1<-name.check(phy,d1)
	td2<-name.check(phy,d2)

	if(td1 != "OK" | td2 != "OK") stop("Names don't match")
	
	# this changes the discrete data to 1:n and remembers the original charStates
	dat1<-as.factor(d1)
	charStates1<-levels(dat1)
	k1<-nlevels(dat1)
	
	dat2<-as.factor(d2)
	charStates2<-levels(dat2)
	k2<-nlevels(dat2)
	
	if(k1 != 2 | k2 != 2) stop("Only 2-state characters currently supported")
	
	mergedDat<-combineDiscreteCharacters(dat1, dat2)
			
	charStates<-levels(mergedDat)
	k<-nlevels(mergedDat)
		
	ndat<-as.numeric(mergedDat)
	names(ndat)<-names(mergedDat)
		
	constraint<-makeDiscreteCorrelationConstraints(modelType="ER")
	
	lik<-make.mkn(phy, ndat, k=k)
	ulik<-constrain(lik, formulae=constraint$uCon, extra=constraint$uExtra)
	clik<-constrain(lik, formulae=constraint$cCon, extra=constraint$cExtra)
	
	unames<-argnames(ulik)
	uML<-find.mle(ulik, setNames(rep(1,length(unames)), argnames(ulik)))
	
	cnames<-argnames(clik)
	cML<-find.mle(clik, setNames(rep(1,length(cnames)), argnames(clik)))
	
	lrStat<-2*(cML$lnLik - uML$lnLik)
	lrDF <- 2*(length(cnames)-length(unames))
	
	lrPVal <- pchisq(lrStat, lrDF, lower.tail=F)
	
	return(list(unames=unames, cnames=cnames, lrStat= lrStat, lrDF= lrDF, lrPVal= lrPVal))

}

combineDiscreteCharacters<-function(dat1, dat2) {
	newDat<-paste(as.character(dat1), as.character(dat2), sep="")
	newDat<-factor(newDat, levels=c("11", "12", "21", "22"))
	names(newDat)<-names(dat1)
	return(newDat)
}
