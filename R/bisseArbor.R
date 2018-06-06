#' BiSSE in aRbor
#'
#' This function tests for a relationship between a two-state character and speciation / extinction rates
#'
#' @param td An object of class 'treedata'
#' @export

bisseArbor<-function(td) {
	
	# check character type
	td <- checkFactor(td, return.factor=TRUE)
	
	res<-list()
	for(i in 1:ncol(td$dat)) {
	  res[[i]] <- bisseArborCalculator(td$phy, td[[i]])
	}
	
	myRes<-list();
	myRes$bisLik<-res[[i]]$fullm$lnLik
	myRes$bisk<-length(res[[i]]$fullm$par)
	myRes$nullLik<-res[[i]]$nullm$lnLik
	myRes$nullk<-length(res[[i]]$nullm$par)
	
	myRes$lrStat <- 2*(myRes$bisLik-myRes$nullLik)
	myRes$lrDF <- myRes$bisk- 	myRes$nullk
	myRes$lrPVal <- pchisq(myRes$lrStat, myRes$lrDF, lower.tail=F)

	paramTable<-c(res[[i]]$nullm$par, res[[i]]$fullm$par)
	names(paramTable)<-c(paste("null_", names(res[[i]]$nullm$par),sep=""),  paste("bis_", names(res[[i]]$fullm$par),sep=""))
	return(list(stats=myRes, param=as.list(paramTable)))
	
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
	
	fm<-getBisseMk2(phy, ndat)
	nm<-getBisseMk2Null(phy, ndat)
	
	zz<-list(fullm=fm, nullm=nm)
	return(zz)	
	
	
}

getBisseMk2<-function(phy, ndat) {
	
	lik<-make.bisse(phy, ndat-1)
	
	p<-starting.point.bisse(phy)
	fit<-find.mle(lik, p)
	
	return(fit)	
}

getBisseMk2Null<-function(phy, ndat) {
  
  lik<-make.bisse(phy, ndat-1)
  
  lik<-constrain(lik, lambda0~lambda1, mu0~mu1)
  
  p<-starting.point.bisse(phy)[argnames(lik)]
  
  fit<-find.mle(lik, p)
  
  return(fit)	
}


			