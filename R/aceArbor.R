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

aceArbor<-function(td, charType="continuous", aceType="marginal", discreteModelType="ER", plot=T) {
	
	# need checks here that the data actually make sense
	if(charType=="continuous") {
		if(!checkNumeric(td)) stop("Data contains factors, which cannot be used for continuous ancestral state reconstruction")
	}
	if(charType=="discrete") {
		if(!checkFactor(td)) warning("Data contains numeric entries, which will be converted to factors for discrete ancestral state reconstruction")
	}
	
	
	res<-lapply(td$dat, function(x) aceArborCalculator(td$phy, x, charType, aceType, discreteModelType, plot))
	
	# the following works for charType = continuous and acetype = marginal
	if(charType=="continuous") {
		ancestralStates<-matrix(nrow=dim(td$dat)[1]-1, ncol=dim(td$dat)[2])
		rownames(ancestralStates)<-names(res[[1]]$ace)
		colnames(ancestralStates)<-colnames(td$dat)
	
		ancestralStatesUpperCI<-ancestralStates
		ancestralStatesLowerCI<-ancestralStates
		
		if(aceType=="MCMC") bayesOutput<-list()

		for(i in 1:length(res)) {
			ancestralStates[,i]<-res[[i]]$ace
			ancestralStatesUpperCI[,i]<-res[[i]]$CI95[,1]
			ancestralStatesLowerCI[,i]<-res[[i]]$CI95[,2]
			if(aceType=="MCMC") bayesOutput[[i]]<-res[[i]]$bayesOutput
		}
	
		if(aceType=="MCMC") {
			res<-list(ancestralStates= ancestralStates, ancestralStatesUpperCI= ancestralStatesUpperCI, ancestralStatesLowerCI= ancestralStatesLowerCI, bayesOutput=bayesOutput)
		} else {
			res<-list(ancestralStates= ancestralStates, ancestralStatesUpperCI= ancestralStatesUpperCI, ancestralStatesLowerCI= ancestralStatesLowerCI)
		}
		return(res)
		
	} else if(charType=="discrete" & acetype=="marginal"){
		
		
	} else return(res)
		
}	

aceArborCalculator<-function(phy, dat, charType="continuous", aceType="marginal", discreteModelType="ER", plot=T, mcmcGen=10000, mcmcBurnin=1000) {
	
	# this function requires a phylo object
 	# and a dat
	# and a colID that tells which column to use
  
	# optional arguments:
	# charType allows the user to force the data to be treated as continuous or discrete.
	# otherwise, factors are treated as discrete and anything else is treated as continuous
	# it might be nice to change this so that continuous data that only has a small number of values is discrete e.g. 0 and 1
	
	# "model" is there for users to pass through details of model specificiation.
	# e.g. OU for continuous, or sym/ard for discrete
	# this is not yet implemented
	
	# check character type
	ctype = match.arg(charType, c("discrete", "continuous"))
	discreteModelType = match.arg(discreteModelType, c("ER", "SYM", "ARD"))
	aceType = match.arg(aceType, c("marginal", "joint", "MCMC"))
	
	if(ctype=="discrete") {
		
		# this changes the discrete data to 1:n and remembers the original charStates
		fdat<-as.factor(dat)
		charStates<-levels(fdat)
		k<-nlevels(fdat)
		
		ndat<-as.numeric(fdat)
		
		if(aceType=="marginal") {
			zz<- getDiscreteAceMarginal(phy, ndat, k, discreteModelType);
		} else if(aceType=="joint") { # this should be modified to average over many reps
			zz<- getDiscreteAceJoint(phy, ndat, k, discreteModelType)
		} else if(aceType=="MCMC"){
			zz<- getDiscreteAceMCMC(phy, ndat, k, discreteModelType)
		}
		
		if(plot) plotDiscreteReconstruction(phy, zz, dat, charStates)
    
		colnames(zz)<-charStates
		return(zz)	
			
	} else if(ctype=="continuous") {
		if(aceType=="marginal") {
			zz<-fastAnc(phy, dat, CI=T)
			names(dat)<-phy$tip.label 
			phenogram(phy, dat)
			return(zz)
		} else if (aceType=="MCMC") {
			names(dat)<-phy$tip.label 
			bayesOutput<-anc.Bayes(phy, dat, ngen= mcmcGen)
			bayesChar<-bayesOutput[,-which(colnames(bayesOutput) %in% c("gen", "sig2", "logLik"))]
			aceStates<-apply(bayesChar, 2, mean)
			CI95 <-t(apply(bayesChar, 2, function(x) quantile(x, c(0.025, 0.975))))
			zz<-list(ace=aceStates, CI95= CI95, bayesOutput=bayesOutput)
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
	


			