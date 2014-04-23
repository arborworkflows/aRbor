library(diversitree)
library(geiger)

aceArbor<-function(phy, dat, charType="fromData", aceType="marginal", discreteModelType="ER") {
	
	# this function requires a tree in ape phylo format (phy)
	# and a single character (dat) with names names(dat) that match phy$tip.label
	
	# optional arguments:
	# charType allows the user to force the data to be treated as continuous or discrete.
	# otherwise, factors are treated as discrete and anything else is treated as continuous
	# it might be nice to change this so that continuous data that only has a small number of values is discrete e.g. 0 and 1
	
	# "model" is there for users to pass through details of model specificiation.
	# e.g. OU for continuous, or sym/ard for discrete
	# this is not yet implemented
	
	# Check that data is a vector
	
	# unlike some functions in geiger, I think we should assume that we have a single character
	# and not a matrix of characters. So I will put in a check. Selecing columns should happen 
	# upstream of this function
	if(!is.null(dim(dat))) stop("This function, aceArbor, requires a vector not a matrix, as input.")
	
	# use treedata to make sure tree and data match
	td<-treedata(phy,dat)
	phy = td$phy
    dat = td$data[,1]
	
	# check character type
	ctype = match.arg(charType, c("fromData", "discrete", "continuous"))
	discreteModelType = match.arg(discreteModelType, c("ER", "SYM", "ASD"))
	aceType = match.arg(aceType, c("marginal", "joint", "MCMC"))
	
	if(ctype=="fromData") # then try to figure it out
	{
		if(is.factor(dat)) {
			charType<-"discrete"
		} else if(nlevels(as.factor(dat))/length(dat) < 0.1) {
			warning("Guessing that this is a discrete character based on repeated values")
			charType<-"discrete"
		} else {
			charType<-"continuous"
		}
	} # needless to say, this is not yet robust
	
	if(ctype=="discrete") {
		dat<-as.factor(dat)
		charStates<-levels(dat)
		k<-nlevels(dat)
		
		ndat<-as.numeric(dat)
		names(ndat)<-names(dat)
		
		lik<-make.mkn(phy, ndat, k=k)
		con<-makeMkConstraints(k=k, modelType= discreteModelType)
		lik<-constrain(lik, q12~q21)
		
		if(aceType=="marginal") {
			fit<-find.mle(lik, setNames(1, argnames(lik)))
			zz<-asr.marginal(lik, coef(fit))
			plot(phy, main="Marginal ASR")
			nodelabels(pie=t(zz), piecol=1:2, cex=.5, frame="circle")
			legend("bottomleft", fill=1:2, legend=charStates)
		} else if(aceType=="joint") { # this should be modified to average over many reps
			fit<-find.mle(lik, setNames(1, argnames(lik)))
			zz<-asr.joint(lik, coef(fit))
			
			plot(phy, main="Joint ASR")
			nodelabels(pch=19, col=zz, bg=zz, cex=2)
			legend("bottomleft", fill=1:2, legend=charStates)
		} else if(aceType=="MCMC"){
			set.seed(1)
			prior <- make.prior.exponential(.5) # NOT GENERAL
			samples <- mcmc(lik, pars, 1000, w=1, prior=prior, print.every=10)
			aceSamp <- apply(samples[c(-1, -dim(samples)[2])], 1, asr.joint, lik=lik)
			zz<-apply(aceSamp, 1, table)/1000

			plot(exTree, main="MCMC ASR")
			nodelabels(pie=t(zz), piecol=1:2, cex=.5, frame="circle")
		}
		
		return(t(zz))
		
	} else if(ctype=="continuous") {
		if(aceType=="marginal") {
			zz<-fastAnc(phy, y, CI=T) 
			phenogram(phy, y)
			return(zz)
		} else if (aceType=="MCMC") {
			zz<-ace.Bayes(phy, y, ngen=10000)
			return(zz)
		} else {
			stop("Not supported yet")
		}
		
	} else stop("Invalid character type in aceArbor.\n")
	
}
