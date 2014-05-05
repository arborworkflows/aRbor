# function for calculating phylogenetic signal of discrete and continuous traits
# can use either Blomberg (for continuous), Pagel lambda (for both), or "garbage test" (for discrete)

physigArbor<-function(phy, dat, charType="fromData", signalTest="pagelLambda", discreteModelType="ER") {
	ctype = match.arg(charType, c("fromData", "discrete", "continuous"))

	if(ctype=="fromData") # then try to figure it out
		ctype<-detectCharacterType(dat)
	
	signalTest = match.arg(signalTest, c("pagelLambda", "Blomberg", "garbageTest"))
	discreteModelType = match.arg(discreteModelType, c("ER", "SYM", "ARD"))

	
	if(ctype=="discrete") {
		# this changes the discrete data to 1:n and remembers the original charStates
		dat<-as.factor(dat)
		charStates<-levels(dat)
		k<-nlevels(dat)
		
		ndat<-as.numeric(dat)
		names(ndat)<-names(dat)
		
		if(signalTest=="pagelLambda") {
			res<-discreteLambdaTest(phy, ndat, k, discreteModelType)
		} else if(signalTest=="garbageTest") {
			
		} else if(signalTest=="Blomberg") {
			stop("Blomberg K should not be used for discrete characters")
		}
	}
	
	if(ctype=="continuous") {
		
	}
	
	return(res)

}

discreteLambdaTest<-function(phy, ndat, k, discreteModelType) {
	m1<-fitDiscrete(phy, ndat, model=discreteModelType)
	m2<-fitDiscrete(phy, ndat, model=discreteModelType, transform="lambda")
	
	chisqTestStat <- 2 * (m2$opt$lnL-m1$opt$lnL)
	chisqPVal <- pchisq(chisqTestStat, 1, lower.tail=F)
	
	aicScores<-c(m1$opt$aicc, m2$opt$aicc)
	names(aicScores)<-c("Mk", "Mk plus lambda")
	
	res<-list(chisqTestStat= chisqTestStat, chisqPVal= chisqPVal, aicScores= aicScores)
	return(res)
}