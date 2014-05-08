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
	
	# key: 1 = 11, 2 = 12, 3 = 21, 4 = 22
	if(modelType=="ER" | modelType == "SYM") {
		uCon<-c("q13~qA", "q31~qA", "q24~qA", "q42~qA", 
				"q12~qB", "q21~qB", "q34~qB", "q43~qB",
				"q14~0", "q41~0", "q23~0", "q32~0")
		uExtra<-c("qA", "qB")
		cCon<-c("q13~qA_B1", "q31~qA_B1", "q24~qA_B2", "q42~qA_B2", 
				"q12~qB_A1", "q21~qB_A1", "q34~qB_A2", "q43~qB_A2",
				"q14~0", "q41~0", "q23~0", "q32~0")
		cExtra<-c("qA_B1","qA_B2", "qB_A1", "qB_A2")

	} else if(modelType=="ARD") {
		uCon<-c("q13~qA12", "q31~qA21", "q24~qA12", "q42~qA21", 
				"q12~qB12", "q21~qB21", "q34~qB12", "q43~qB21",
				"q14~0", "q41~0", "q23~0", "q32~0")
		uExtra<-c("qA12", "qA21", "qB12", "qB21")
		cCon<-c("q13~qA12_B1", "q31~qA21_B1", "q24~qA12_B2", "q42~qA21_B2", 
				"q12~qB12_A1", "q21~qB21_A1", "q34~qB12_A2", "q43~qB21_A2",
				"q14~0", "q41~0", "q23~0", "q32~0")
		cExtra<-c("qA12_B1","qA21_B1","qA12_B2","qA21_B2", "qB12_A1","qB21_A1", "qB12_A2", "qB21_A2")
		
	}
	
	return(list(uCon=uCon, uExtra=uExtra, cCon=cCon, cExtra=cExtra))
	
}
