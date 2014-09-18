#' Carry out a phylogenetic PCA analysis
#'
#' @param td An object of class 'treedata'
#' @export


phylopcaArbor<-function(td) {
	
  td <- checkNumeric(td, return.numeric=TRUE)
  ncol<-dim(td$dat)[2]

  res<-phyl.pca(td$phy, td$dat)

  return(res)

}