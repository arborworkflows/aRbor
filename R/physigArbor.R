#' Function for calculating phylogenetic signal of discrete and continuous traits
#' 
#' This function allows testing either Blomberg (for continuous), Pagel lambda (for both), or "garbage test" (for discrete)
#' @export

physigArbor<-function(td, charType="fromData", signalTest="pagelLambda", discreteModelType="ER") {
	charType = match.arg(charType, c("fromData", "discrete", "continuous"))

	if(charType =="fromData") # then try to figure it out
		charType <-detectCharacterType(td$dat)
	
	# check that the data actually make sense - this is a pretty weak test
	if(charType=="continuous") {
    	td <- checkNumeric(td, return.numeric=TRUE)
	}
	if(charType=="discrete") {
		td <- checkFactor(td, return.factor=TRUE)
	}
	
	signalTest = match.arg(signalTest, c("pagelLambda", "Blomberg", "garbageTest"))
	discreteModelType = match.arg(discreteModelType, c("ER", "SYM", "ARD"))

	if(any(is.na(td$dat))){
    	if(na.rm=="bytrait"){
      		res <- lapply(1:ncol(td$dat), function(i) {
        	tdi <- select(td, i);
       		tdi <- filter(tdi, !is.na(tdi$dat[,1]));
        	physigArborCalculator(tdi$phy, setNames(tdi$dat[,1], rownames(tdi$dat)), charType, signalTest, discreteModelType)
      		})
    	}
    	if(na.rm=="all"){
      		td <- filter(td, apply(apply(td$dat, 1, function(x) !is.na(x)), 2, all))
      		res <- lapply(td$dat, function(x) physigArborCalculator(td$phy, setNames(x, rownames(td$dat)), charType, signalTest, discreteModelType))
    	}
  	} else {
	  res <- lapply(td$dat, function(x) physigArborCalculator(td$phy, setNames(x, rownames(td$dat)), charType, signalTest, discreteModelType))
  	}
  	
	class(res) <- c("physigArbor", class(res))
  	attributes(res)$td <- td
	attributes(res)$charType <- charType
	attributes(res)$signalTest <- signalTest
  	if(charType=="discrete"){
    	attributes(res)$discreteModelType = discreteModelType
    	attributes(res)$charStates = lapply(1:ncol(td$dat), function(x) levels(td$dat[,x]))

  	}
  	if(any(is.na(td$dat))){
    	attributes(res)$na.drop <- lapply(td$dat, function(x) rownames(td$dat)[which(is.na(x))])
  	}
	names(res) <- colnames(td$dat)
	return(res)
}

physigArborCalculator<-function(phy, dat, charType="fromData", signalTest="pagelLambda", discreteModelType="ER") {
	
	if(charType =="discrete") {
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
	
	if(charType =="continuous") {
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
	
	lnlValues<-c(m1$opt$lnL, m2$opt$lnL)
	names(lnlValues)<-c("Lambda fixed at one", "Lambda estimated")
	
	aiccScores<-c(m1$opt$aicc, m2$opt$aicc)
	names(aiccScores)<-c("Lambda fixed at one", "Lambda estimated")
	
	res<-list(lnlValues= lnlValues, chisqTestStat= chisqTestStat, chisqPVal= chisqPVal, aiccScores= aiccScores)
	return(res)
}

discreteGarbageTest<-function(phy, ndat, k, discreteModelType) {
	m1<-fitDiscrete(phy, ndat, model=discreteModelType)
	m2<-fitDiscreteGarbageModel(phy, ndat)
	
	lnlValues<-c(m1$opt$lnL, m2$lnL)
	names(lnlValues)<-c("Mk", "Garbage")
	
	
	aiccScores<-c(m1$opt$aicc, m2$aicc)
	names(aiccScores)<-c("Mk", "Garbage")
	
	res<-list(lnlValues= lnlValues, chisqTestStat= NULL, chisqPVal= NULL, aiccScores= aiccScores)
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
	
	lnlValues<-c(m1$opt$lnL, m2$opt$lnL)
	names(lnlValues)<-c("BM", "BM+lambda")
	
	
	chisqTestStat <- 2 * (m2$opt$lnL-m1$opt$lnL)
	chisqPVal <- pchisq(chisqTestStat, 1, lower.tail=F)
	
	aiccScores<-c(m1$opt$aicc, m2$opt$aicc)
	names(aiccScores)<-c("BM", "BM+lambda")
	
	res<-list(lnlValues = lnlValues, chisqTestStat= chisqTestStat, chisqPVal= chisqPVal, aiccScores= aiccScores)
	return(res)
}


continuousGarbageTest<-function(phy, dat) {
	m1<-fitContinuous(phy, dat)
	m2<-fitContinuous(phy, dat, model="white")
	
	lnlValues<-c(m1$opt$lnL, m2$opt$lnL)
	names(lnlValues)<-c("BM", "white-noise")
	
	
	
	aiccScores<-c(m1$opt$aicc, m2$opt$aicc)
	names(aiccScores)<-c("Mk", "WhiteNoise")
	
	res<-list(lnlValues= lnlValues, chisqTestStat= NULL, chisqPVal= NULL, aiccScores= aiccScores)
	return(res)
}

continuousBlombergTest<-function(phy, dat) {
	pstest<-phylosignal(dat, phy)
	bres<-list(k= pstest$K, vObs= pstest$PIC.variance.obs, vRnd=pstest$PIC.variance.rnd.mean, pVal=pstest$PIC.variance.P, zScore=pstest$PIC.variance.Z)
	res<-list(lnlValues= NULL, chisqTestStat= NULL, chisqPVal= NULL, aiccScores= NULL, bres=bres)
	return(res)
}

#' @export
print.physigArbor <- function(x, ...){
  names <- attributes(x)$names
  attributes(x) <- NULL
  attributes(x)$names <-  names
  print(x)
}
