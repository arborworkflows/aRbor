#' Function for calculating phylogenetic signal of discrete and continuous traits
#' 
#' This function allows testing either Blomberg (for continuous), Pagel lambda (for both), or "garbage test" (for discrete)
#' @export

physigArbor<-function(phy, dat, charType="fromData", signalTest="pagelLambda", discreteModelType="ER") {
	ctype = match.arg(charType, c("fromData", "discrete", "continuous"))

	if(ctype=="fromData") # then try to figure it out
		ctype<-detectCharacterType(dat)
	
	signalTest = match.arg(signalTest, c("pagelLambda", "Blomberg", "garbageTest"))
	discreteModelType = match.arg(discreteModelType, c("ER", "SYM", "ARD"))

	missingSpecies<-is.na(dat)
	if(sum(missingSpecies) > 0) {
		dat<-dat[!missingSpecies]
		phy<-drop.tip(phy, names(dat)[missingSpecies])
	}
	
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
			res<-discreteGarbageTest(phy, ndat, k, discreteModelType)
		} else if(signalTest=="Blomberg") {
			stop("Blomberg K should not be used for discrete characters")
		}
	}
	
	if(ctype=="continuous") {
		if(signalTest=="pagelLambda") {
			res<-continuousLambdaTest(phy, dat)
		} else if(signalTest=="garbageTest") {
			res<-continuousGarbageTest(phy, dat)
		} else if(signalTest=="Blomberg") {
			res<-continuousBlombergTest(phy, dat)
		}
	}
	
	return(res)

}

discreteLambdaTest<-function(phy, ndat, k, discreteModelType) {
	m1<-fitDiscrete(phy, ndat, model=discreteModelType)
	m2<-fitDiscrete(phy, ndat, model=discreteModelType, transform="lambda")
	
	chisqTestStat <- 2 * (m2$opt$lnL-m1$opt$lnL)
	chisqPVal <- pchisq(chisqTestStat, 1, lower.tail=F)
	
	aicScores<-c(m1$opt$aicc, m2$opt$aicc)
	names(aicScores)<-c("Mk", "Mk+lambda")
	
	res<-list(chisqTestStat= chisqTestStat, chisqPVal= chisqPVal, aicScores= aicScores)
	return(res)
}

discreteGarbageTest<-function(phy, ndat, k, discreteModelType) {
	m1<-fitDiscrete(phy, ndat, model=discreteModelType)
	m2<-fitDiscreteGarbageModel(phy, ndat)
	
	aiccScores<-c(m1$opt$aicc, m2$aicc)
	names(aiccScores)<-c("Mk", "Garbage")
	
	res<-list(aiccScores= aiccScores)
	return(res)
}

fitDiscreteGarbageModel<-function(phy, ndat) {
	tt <- table(ndat)
	prob <- tt/sum(tt)
	lnL <- sum(tt * log(prob))
	k <- length(tt)-1
	n <- length(ndat)
	aic <- -2 * lnL + 2 * k
	aicc <- aic + 2 * k * (k+1) / (n - k - 1)
	res<-list(prob=prob, lnL=lnL, k=k, aic=aic, aicc=aicc)
	return(res)
}


continuousLambdaTest<-function(phy, dat) {
	m1<-fitContinuous(phy, dat, model="BM")
	m2<-fitContinuous(phy, dat, model="lambda")
	
	chisqTestStat <- 2 * (m2$opt$lnL-m1$opt$lnL)
	chisqPVal <- pchisq(chisqTestStat, 1, lower.tail=F)
	
	aicScores<-c(m1$opt$aicc, m2$opt$aicc)
	names(aicScores)<-c("BM", "BM+lambda")
	
	res<-list(chisqTestStat= chisqTestStat, chisqPVal= chisqPVal, aicScores= aicScores)
	return(res)
}


continuousGarbageTest<-function(phy, dat) {
	m1<-fitContinuous(phy, dat)
	m2<-fitContinuous(phy, dat, model="white")
	
	aiccScores<-c(m1$opt$aicc, m2$opt$aicc)
	names(aiccScores)<-c("Mk", "WhiteNoise")
	
	res<-list(aiccScores= aiccScores)
	return(res)
}

continuousBlombergTest<-function(phy, dat) {
	res<-phylosignal(dat, phy)
	return(res)
}