# function for calculating phylogenetic signal of discrete and continuous traits
# can use either Blomberg (for continuous), Pagel lambda (for both), or "garbage test" (for discrete)

aceArbor<-function(phy, dat, charType="fromData", signalTest="pagelLambda", discreteModelType="ER") {
	ctype = match.arg(charType, c("fromData", "discrete", "continuous"))

	if(ctype=="fromData") # then try to figure it out
		ctype<-detectCharacterType(dat)
	
	signalTest = match.arg(aceType, c("pagelLambda", "Blomberg", "garbageTest"))

	
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

}

discreteLambdaTest<-function(phy, ndat, k, discreteModelType) {
	m1<-fitDiscrete(phy, ndat, model=discreteModelType)
	m2<-fitDiscrete(phy, ndat, model= discreteModelType, transform="lambda")
}