#' BiSSE in aRbor
#'
#' This function tests for a relationship between a two-state character and speciation / extinction rates
#'
#' @param td An object of class 'treedata'
#' @export

bisseArbor<-function(td) {
	
	# check character type
	td <- checkFactor(td, return.factor=TRUE)
	
	res <- lapply(td$dat, function(x) bisseArborCalculator(td$phy, setNames(x, rownames(td$dat))))
	
	class(res) <- c("asrArbor", class(res))
  	attributes(res)$td <- td
  	attributes(res)$charType <- "discrete"
	attributes(res)$bisseType <- "Mk2"
    attributes(res)$charStates = lapply(1:ncol(td$dat), function(x) levels(td$dat[,x]))
 
	return(res)
	
}	

bisseArborCalculator<-function(phy, dat, names=NULL) {
	
	# this function requires a phylo object
 	# and a dat
		
	# this changes the discrete data to 1:n and remembers the original charStates
	fdat<-as.factor(dat)
	charStates<-levels(fdat)
	#k<-nlevels(fdat)
		
	ndat<-as.numeric(fdat)
    names(ndat) <- names(fdat)
		
	if(!is.null(names)) names(ndat)<-names
	
	zz<-getBisseMk2(phy, ndat)
	
	return(zz)	
	
	
}

getBisseMk2<-function(phy, ndat) {
	
	lik<-make.bisse(phy, ndat-1)
	
	#model type not supported yet
	#con<-makeMkConstraints(k=k, modelType= discreteModelType)
	#ltemp<-lik
	
	#if(!is.null(con))
	#	for(i in 1:length(con)) ltemp<-constrain(ltemp, con[[i]], extra=extra)
	
	#lik<-ltemp
	
	p<-starting.point.bisse(phy)
	fit<-find.mle(lik, p)
	
	return(fit)	
}



			