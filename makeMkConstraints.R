makeMkConstraints<-function(k, modelType="ER") {
	  modelType = match.arg(modelType, c("ER", "SYM", "ARD", "other"))
	  
	  res<-character()
	  
	  if(modelType=="SYM") {
	  	for(i in 1:(k-1))
	  		for(j in (i+1):k)
	  			res<-c(res, paste("q",i,j,"~q",j,i, sep=""))
	  			
	  } else if(modelType=="ER") {
	  	for(i in 1:(k-1))
	  		for(j in (i+1):k)
	  			res<-c(res, paste("q",i,j,"~q",j,i, sep=""))
	  	
	  	for(j in 3:k)
	  		res<-c(res, paste("q",1,2,"~q",1,j, sep=""))
	  				
	  } else if(modelType=="ARD") {
	  	
	  } else if(modelType=="other") {
	  	stop("Not supported")
	  } else {
	  	stop("Error specifying Mk model type")
	  }
	res
}