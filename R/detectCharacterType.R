detectCharacterType<-function(dat, cutoff=0.1) {
	if(is.factor(dat)) {
			charType<-"discrete"
	} else if(nlevels(as.factor(dat))/length(dat) < cutoff) {
			warning("Guessing that this is a discrete character based on repeated values")
			charType<-"discrete"
	} else {
			charType<-"continuous"
	}
	return(charType)
}		# needless to say, this is not yet robust