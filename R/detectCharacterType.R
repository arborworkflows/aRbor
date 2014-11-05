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

# Row and column name check and force
hasNames <- function(dat, nameType="row") {
	nType = match.arg(nameType, c("row", "col", "rowcol"))
	if(nType=="row") {
		res<-!is.null(rownames(dat))
	}
	if(nType=="col") {
		res<-!is.null(colnames(dat))
	}
	if(nType=="rowcol") {
		res<-!is.null(rownames(dat)) & !is.null(colnames(dat))
	}
	names(res)<-nameType
	res
}

forceNames <- function(dat, nameType="row") {
	nType = match.arg(nameType, c("row", "col", "rowcol"))
	if(nType=="row" | nType=="rowcol") {
		if(!hasNames(dat, nameType="row")) {
			nrows<-dim(dat)[1]
			rownames(dat) <- paste("n", 1:nrows, sep="")
		}	
	}
	if(nType=="col" | nType=="rowcol") {
		if(!hasNames(dat, nameType="col")) {
			ncols<-dim(dat)[2]
			colnames(dat) <- paste("n", 1:ncols, sep="")
		}		
	}

	dat
}
