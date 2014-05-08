#' Reconstructing ancestral states of discrete and continuous characters
#'
#' This function returns ancestral state estimates using a variety of methods
#'
#' @param td An object of class 'treedata'
#' @param colID A column selector for the dataframe in td
#' @param charType specifies the type of character, either:
#' \describe{
#' 		\item{"discrete"}{a character with a discrete number of states}
#' 		\item{"continuous"}{a continuously varying character}	
#' 		\item{"fromData"}{will attempt to determine the data type from the data itself; DO NOT USE YET!}
#'	}
#' @param aceType specifies the method used to reconstruct ancestral character states:
#' \describe{
#' 		\item{"marginal"}{marginal ancestral state reconstructions, which reconstruct each node integrating over all possibilities at all other nodes in the tree; this is typically the method used in the literature to reconstruce ACEs}
#' 		\item{"joint"}{joint ancestral reconstructions, which give the configuration of ancestral states that together maximize the likelihood of the data given model parameters}	
#' 		\item{"mcmc"}{reconstruct ancestral states using Bayesian MCMC. Note that the discrete version of this doesn't seem to work, and even if it did work it is not a full MCMC ancestral state method}
#'	}	
#' @param discreteModelType One of ER, SYM, or ARD; see geiger's fitDiscrete for full description
#' @param plot If true, make a plot of ancestral states.
	

aceArbor<-function(td, colID, charType="fromData", aceType="marginal", discreteModelType="ER", plot=T) {
	
	# this function requires a treedata object
 	# (see make.treedata)
	# and a colID that tells which column to use
  
	# optional arguments:
	# charType allows the user to force the data to be treated as continuous or discrete.
	# otherwise, factors are treated as discrete and anything else is treated as continuous
	# it might be nice to change this so that continuous data that only has a small number of values is discrete e.g. 0 and 1
	
	# "model" is there for users to pass through details of model specificiation.
	# e.g. OU for continuous, or sym/ard for discrete
	# this is not yet implemented
	
  # obtain the column that you want
  td<-select(td, colID)
	
	# check character type
	ctype = match.arg(charType, c("fromData", "discrete", "continuous"))
	discreteModelType = match.arg(discreteModelType, c("ER", "SYM", "ARD"))
	aceType = match.arg(aceType, c("marginal", "joint", "MCMC"))
	
	if(ctype=="fromData") # then try to figure it out
		ctype<-detectCharacterType(td$dat)
	
	if(ctype=="discrete") {
		
		# this changes the discrete data to 1:n and remembers the original charStates
		dat<-as.factor(td$dat[,1])
    names(dat)<-rownames(td$dat)
		charStates<-levels(dat)
		k<-nlevels(dat)
		
		ndat<-as.numeric(dat)
		names(ndat)<-names(dat)
		
		if(aceType=="marginal") {
			zz<- getDiscreteAceMarginal(td$phy, ndat, k, discreteModelType);
		} else if(aceType=="joint") { # this should be modified to average over many reps
			zz<- getDiscreteAceJoint(td$phy, ndat, k, discreteModelType)
		} else if(aceType=="MCMC"){
			zz<- getDiscreteAceMCMC(td$phy, ndat, k, discreteModelType)
		}
		
		if(plot) plotDiscreteReconstruction(td$phy, zz, dat, charStates)
    
		colnames(zz)<-charStates
		return(zz)	
			
	} else if(ctype=="continuous") {
		if(aceType=="marginal") {
			zz<-fastAnc(td$phy, td$dat, CI=T) 
			phenogram(td$phy, td$dat)
			return(zz)
		} else if (aceType=="MCMC") {
			zz<-anc.Bayes(td$phy, td$dat, ngen=10000)
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

plotDiscreteReconstruction<-function(phy, zz, dat, charStates) {
	plot(phy, main="ASR")
	nodelabels(pie=zz, piecol=1:2, cex=.5, frame="circle")
	tiplabels(pch=21, bg=as.numeric(factor(dat))) 
	legend("bottomleft", fill=1:2, legend=charStates)
}


getDiscreteAceJoint<-function(phy, ndat, k, discreteModelType) {
	reps <- 1000
  lik<-make.mkn(phy, ndat, k=k)
	con<-makeMkConstraints(k=k, modelType= discreteModelType)
	if(!is.null(con))
		lik<-constrain(lik, con)
				
	pnames<-argnames(lik)
	fit<-find.mle(lik, setNames(rep(1,length(pnames)), argnames(lik)))
	
	xx <-sapply(1:reps, function(x) asr.joint(lik, coef(fit)))
  zz<-matrix(0, nrow=length(xx), ncol=k)
	zz[,1]<- apply(xx, 1, function(x) sum(x==1)/reps)
	zz[,2]<- apply(xx, 1, function(x) sum(x==2)/reps)
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
	


			