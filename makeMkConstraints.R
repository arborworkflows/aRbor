makeMkConstraints<-function(k, modelType="ER") {
	  modelType = match.arg(modelType, c("ER", "SYM", "ARD", "other"))
	  
	  res<-character()
	  
	  if(modelType=="SYM" | modelType=="ER") {
	  	for(i in 1:(k-1))
	  		for(j in (i+1):k)
	  			res<-c(res, paste("q",i,j,"~q",j,i, sep=""))
	  			
	  } 
	  
	  if(modelType=="ER") {
	  	for(j in 3:k)
	  		for(i in 1:(j-1))
	  			res<-c(res, paste("q",1,2,"~q",i,j, sep=""))
	  		  				
	  } 
	  
	  if(modelType=="ARD") {
	  	res<-NULL
	  } 
	  
	  if(modelType=="other") {
	  	stop("Not supported")
	  } 
	  
	 
	res
}

# this function is not general at all - it only works for k=2 and ER model
makeDiscreteCorrelationConstraints<-function(modelType="ER") {
	
	# correlated model - both depend on one another
	# cCon<-makeMkConstraints(k=4, model= modelType)
	# key: 1 = 11, 2 = 12, 3 = 21, 4 = 22
	if(modelType=="ER" | modelType == "SYM") {
		uCon<-c("q13~qC1", "q31~qC1", "q24~qC1", "q42~qC1", 
				"q12~qC2", "q21~qC2", "q34~qC2", "q43~qC2",
				"q14~0", "q41~0", "q23~0", "q32~0")
		uExtra<-c("qC1", "qC2")
		cCon<-c("q13~qC11", "q31~qC11", "q24~qC12", "q42~qC12", 
				"q12~qC21", "q21~qC21", "q34~qC22", "q43~qC22",
				"q14~0", "q41~0", "q23~0", "q32~0")
		cExtra<-c("qC11","qC12", "qC21", "qC22")

	}
	# Add ARD model
	
	return(list(uCon=uCon, uExtra=uExtra, cCon=cCon, cExtra=cExtra))
	
}
