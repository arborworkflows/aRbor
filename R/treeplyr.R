#' Function for making an object of class 'treedata'
#' 
#' This function generates an object of class 'treedata' that ensures that the ordering of tip labels and data remain intact. The object can be manipulated using \code{dplyr} functions.
#' 
#' @param tree An object of class 'phylo'
#' @param data A data frame or matrix
#' @param name_column The column of \code{data} that contains the names to be matched to the tree. By default, it is set to 0, which indicates row names.
#' @details
#' @return An object of class "\code{treedata}". The tree is pruned of tips not represented in the data, and the data is filtered for taxa not in the tree. The data is returned as a data frame tble that is compatible with \code{dplyr} functions. 
#' @export
make.treedata <- function(tree, data, name_column=0) {
  dat <- tbl_df(data)
  if(name_column==0){
    clnm <- colnames(dat)
    dat <- mutate(dat, tip.label=rownames(dat))
    dat <- dat[,c("tip.label", clnm)]
  } else {
    clnm <- ifelse(is.numeric(name_column), 1:ncol(dat)[-name_column], colnames(dat)[-which(colnames(dat)==name_column)])
    dat <- dat[, c(name_column, clnm)]
  }
  data_not_tree <- setdiff(dat[,1], tree$tip.label)
  tree_not_data <- setdiff(tree$tip.label, dat[,1])
  phy <- drop.tip(tree, tree_not_data)
  dat <- filter(dat, tip.label %in% phy$tip.label)
  o <- match(dat[,1], phy$tip.label)
  dat <- arrange(dat, o)
  td <- list(phy=phy, data=dat)
  class(td) <- c("treedata", "list")
  return(td)
}

#' @export
mutate.treedata <- function(tdObject, ...){
  dat <- mutate(tdObject$data, ...)
  tdObject$data <- dat
  return(tdObject)
}

#' @export
select.treedata <- function(tdObject, ...){
  dat <- select(tdObject$data, ...)
  #dat <- lapply(dat, function(x){names(x) <- tdObject$data[,1]; x})
  rownames(dat) <- tdObject$data[,1]
  tdObject$data <- dat
  return(tdObject)
}

#' @export
filter.treedata <- function(tdObject, ...){
  tdObject$data <- filter(tdObject$data, ...)
  tdObject$phy <- drop.tip(tdObject$phy, tdObject$phy$tip.label[!(tdObject$phy$tip.label %in% tdObject$data[,1])])
  return(tdObject)
}

#' @export
summarize.treedata <- function(tdObject, ...){
  res <- summarize(tdObject$data, ...)
  return(res)
}

#' @export
summarise.treedata <- function(tdObject, ...){
  res <- summarise(tdObject$data, ...)
  return(res)
}

