#' Test for an evolutionary correlation between two continuous characters
#'
#' @param td An object of class 'treedata'
#' @param colID A column selector for the dataframe in td


continuousCorrelation<-function(td) {
	
  td <- checkNumeric(td, return.numeric=TRUE)
  ncol<-dim(td$dat)[2]

  vv<-colnames(td$dat)
    
  res<-list()
  counter<-1
  for(i in 1:(ncol-1))
    for(j in (i+1):ncol) {

      xvar<-vv[i]
      yvar<-vv[j]
      
      res[[counter]]<-gls(reformulate(termlabels=xvar, response=yvar), correlation=corBrownian(phy=td$phy), data=td$dat)
      names(res)[counter]<-deparse(reformulate(termlabels=xvar, response=yvar))
      counter<-counter+1
    }
  
	
	return(res)

}