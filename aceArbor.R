library(diversitree)
library(geiger)

aceArbor<-function(phy, dat, charType="fromData", aceType="marginal", discreteModelType="ER", plot=T) {
	
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
	discreteModelType = match.arg(discreteModelType, c("ER", "SYM", "ARD"))
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
		
		# this changes the discrete data to 1:n and remembers the original charStates
		dat<-as.factor(dat)
		charStates<-levels(dat)
		k<-nlevels(dat)
		
		ndat<-as.numeric(dat)
		names(ndat)<-names(dat)
		
		if(aceType=="marginal") {
			zz<-getDiscreteAceMarginal(phy, ndat, k, discreteModelType);
			if(plot) plotDiscreteReconstruction(phy, zz, charStates)
		} else if(aceType=="joint") { # this should be modified to average over many reps
			zz<-getDiscreteAceJoint(phy, ndat, k, discreteModelType)
			if(plot) plotDiscreteReconstruction(phy, zz, charStates)
		} else if(aceType=="MCMC"){
			zz<- getDiscreteAceMCMC(phy, ndat, k, discreteModelType)
			if(plot) plotDiscreteReconstruction(phy, zz, charStates)
		}
		
		colnames(zz)<-charStates
		return(zz)	
			
	} else if(ctype=="continuous") {
		if(aceType=="marginal") {
			zz<-fastAnc(phy, y, CI=T) 
			phenogram(phy, y)
			return(zz)
		} else if (aceType=="MCMC") {
			zz<-anc.Bayes(phy, y, ngen=10000)
			return(zz)
		} else {
			stop("Not supported yet")
		}
		
	} else stop("Invalid character type in aceArbor.\n")
	
}

getDiscreteAceMarginal<-function(phy, ndat, k, discreteModelType) {
	lik<-make.mkn(phy, ndat, k=k)
	con<-makeMkConstraints(k=k, modelType= discreteModelType)
	if(!is.null(con))
		lik<-constrain(lik, con)
	
	pnames<-argnames(lik)
	fit<-find.mle(lik, setNames(rep(1,length(pnames)), argnames(lik)))
	
	zz<-t(asr.marginal(lik, coef(fit)))
	zz		
}

plotDiscreteReconstruction<-function(phy, zz, charStates) {
	plot(phy, main="ASR")
	nodelabels(pie=zz, piecol=1:2, cex=.5, frame="circle")
	legend("bottomleft", fill=1:2, legend=charStates)
}


getDiscreteAceJoint<-function(phy, ndat, k, discreteModelType) {
	lik<-make.mkn(phy, ndat, k=k)
	con<-makeMkConstraints(k=k, modelType= discreteModelType)
	if(!is.null(con))
		lik<-constrain(lik, con)
				
	pnames<-argnames(lik)
	fit<-find.mle(lik, setNames(rep(1,length(pnames)), argnames(lik)))
	
	xx<-asr.joint(lik, coef(fit))
	zz<-matrix(0, nrow=length(xx), ncol=k)
	zz[which(xx==1),1]<-1	
	zz[which(xx==2),2]<-1		
	zz
}

getDiscreteAceMCMC<-function(phy, ndat, k, discreteModelType) { # results do not look correct to me
	lik<-make.mkn(phy, ndat, k=k)
	con<-makeMkConstraints(k=k, modelType= discreteModelType)
	if(!is.null(con))
		lik<-constrain(lik, con)
				
	set.seed(1) # this is not general
	prior <- make.prior.exponential(.5) # NOT GENERAL - we need arguments for this
	pars<-exp(prior(1))
	samples <- mcmc(lik, pars, 1000, w=1, prior=prior, print.every=10) # likewise, need control arguments here
	aceSamp <- apply(samples[c(-1, -dim(samples)[2])], 1, asr.joint, lik=lik)
	zz<-apply(aceSamp, 2, table)/1000
	t(zz)
}
	


			