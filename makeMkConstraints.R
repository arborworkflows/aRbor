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