
#' Test for an evolutionary correlation between two continuous characters
#'
#' @param tree A phylogenetic tree
#' @param table A data matrix
#' @param column1 One column from table; must be a binary discrete trait
#' @param column2 Second column from table; must be a binary discrete trait
#' @param modelType Discrete model, one of "ER", "SYM", "ARD"
#' @export
discreteCorrelation<-function(tree, table, column1, column2, modelType="ER") {

  # check for perfect name matching
  # this could be improved a lot

  td <- make.treedata(tree, table)

  td1 <- select(td, column1, column2)

  # this changes the discrete data to 1:n and remembers the original charStates
  datA<-as.factor(td1[[column1]])
  charStatesA<-levels(datA)
  kA<-nlevels(datA)

  datB<-as.factor(td1[[column2]])
  charStatesB<-levels(datB)
  kB<-nlevels(datB)

  if(kA != 2 | kB != 2) stop("Only 2-state characters currently supported")

  mergedDat<-combineDiscreteCharacters(datA, datB)

  charStates<-levels(mergedDat)
  k<-nlevels(mergedDat)

  ndat<-as.numeric(mergedDat)
  names(ndat)<-names(mergedDat)

  constraint<-makeDiscreteCorrelationConstraints(modelType=modelType)

  lik<-make.mkn(tree, ndat, k=k)
  ulik<-constrain(lik, formulae=constraint$uCon, extra=constraint$uExtra)
  clik<-constrain(lik, formulae=constraint$cCon, extra=constraint$cExtra)

  unames<-argnames(ulik)
  uML<-find.mle(ulik, setNames(rep(1,length(unames)), argnames(ulik)))

  cnames<-argnames(clik)
  cML<-find.mle(clik, setNames(rep(1,length(cnames)), argnames(clik)))

  lrStat<-2*(cML$lnLik - uML$lnLik)
  lrDF <- 2*(length(cnames)-length(unames))

  lrPVal <- pchisq(lrStat, lrDF, lower.tail=F)

  res<-c(lrStat,lrDF,lrPVal, uML$par, cML$par)
  res<-data.frame(c("lrStat","lrDF","lrPVal", unames, cnames), res)
  colnames(res)<-c("Name", "value")
  return(res)

}



combineDiscreteCharacters<-function(dat1, dat2) {
	newDat<-interaction(dat1, dat2, lex.order=T) # lex.order to follow notation order in Pagel
	names(newDat)<-names(dat1)
	return(newDat)
}
