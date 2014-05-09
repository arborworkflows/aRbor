#' Function for making an object of class 'treedata'
#' 
#' This function generates an object of class 'treedata' that ensures that the ordering of tip labels and data remain intact. The object can be manipulated using \code{dplyr} functions.
#' 
#' @param tree An object of class 'phylo'
#' @param data A data frame or matrix
#' @param name_column The column of \code{data} that contains the names to be matched to the tree. By default, it is set to 0, which indicates row names.
#' @return An object of class "\code{treedata}". The tree is pruned of tips not represented in the data, and the data is filtered for taxa not in the tree. The data is returned as a data frame tble that is compatible with \code{dplyr} functions. 
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat, name_column=1)
#' @export
make.treedata <- function(tree, data, name_column=0) {
  coln <- colnames(data)
  if(name_column==0){
      dat.label <- rownames(data)
    } 
  dat <- tbl_df(as.data.frame(lapply(1:ncol(data), function(x) type.convert(as.character(data[,x])))))
  colnames(dat) <- coln
  #dat <- apply(dat, 2, type.convert)
  if(name_column==0){
    clnm <- colnames(dat)
    dat <- dat[,cln]
  } else {
    if(is.numeric(name_column)){
      clnm <- (1:ncol(data))[-name_column]
    } else {
      clnm <- colnames(dat)[-which(colnames(dat)==name_column)]
    }
    dat <- dat[, clnm]
    dat.label <- as.character(data[,name_column])
  }
  data_not_tree <- setdiff(dat.label, tree$tip.label)
  tree_not_data <- setdiff(tree$tip.label, dat.label)
  phy <- drop.tip(tree, tree_not_data)
  dat <- filter(dat, dat.label %in% phy$tip.label)
  o <- match(dat.label, phy$tip.label)
  dat <- arrange(dat, o)
  td <- list(phy=phy, data=dat)
  class(td) <- c("treedata", "list")
  attributes(td)$tip.label <- phy$tip.label
  rownames(td$data) <- attributes(td)$tip.label
  return(td)
}

#' @export
mutate.treedata <- function(tdObject, ...){
  dat <- mutate(tdObject$data, ...)
  tdObject$data <- dat
  rownames(tdObject$data) <- attributes(tdObject)$tip.label
  return(tdObject)
}

#' @export
select.treedata <- function(tdObject, ...){
  dat <- select(tdObject$data, ...)
  #dat <- lapply(dat, function(x){names(x) <- tdObject$data[,1]; x})
  tdObject$data <- dat
  rownames(tdObject$data) <- attributes(tdObject)$tip.label
  return(tdObject)
}


#' @export
filter.treedata <- function(tdObject, ...){
  tdObject$data <- mutate(tdObject$data, tip.label=attributes(tdObject)$tip.label)
  tdObject$data <- filter(tdObject$data, ...)
  attributes(tdObject)$tip.label <- tdObject$data$tip.label
  tdObject$data <- select(tdObject$data, 1:(ncol(tdObject$data)-1))
  rownames(tdObject$data) <- attributes(tdObject)$tip.label
  tdObject$phy <- drop.tip(tdObject$phy, tdObject$phy$tip.label[!(tdObject$phy$tip.label %in% attributes(tdObject)$tip.label)])
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

#' @export
reorder.treedata <- function(tdObject, order="postorder", index.only=FALSE, ...){
  phy <- reorder(tdObject$phy, order, index.only)#, ...)
  index <- match(tdObject$phy$tip.label, phy$tip.label)
  tdObject$dat <- tdObject$dat[index,]
  attributes(tdObject)$tip.label <- phy$tip.label
  tdObject$phy <- phy
  if(index.only){
    return(index)
  } else {
    return(tdObject)
  }
}

#' @export
treeply <- function(tdObject, ...){
  UseMethod("treeply")
}

#' @export
treeply.treedata <- function(tdObject, FUN, ...){
  if(!class(tdObject)[1]=="treedata") stop("Object is not of class 'treedata'")
  FUN <- match.fun(FUN)
  out_FUN <- FUN(tdObject$phy, ...)
  if(class(out_FUN)[1] == "phylo"){
    tdObject$phy <- out_FUN
    tdObject$dat <- tdObject$dat[tdObject$phy$tip.label,]
    attributes(tdObject)$tip.label <- tdObject$phy$tip.label
    return(tdObject)
  } else {
    warning("Function output was not of class 'phylo', returning output only")
    return(out_FUN)
  } 
}

#' @export
treedply <- function(tdObject, ...){
  UseMethod("treedply")
}

#' Run a function on a 'treedata' object
#' 
#' @param tdObject A treedata object
#' @param ... A function call.
#' 
#' @details This function allows arbitrary R functions that use trees and data to be run on a 'treedata' objects. 
#' 
#' @return Function output
#' 
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat, name_column=1)
#' treedply(td, fitContinuous(phy, getVector(dat, SVL), model="BM"))
#' treedply(td, fitContinuous(phy, select(dat, PCI_limbs, PCII_head, 
#'                PCIII_padwidth_vs_tail, PCIV_lamella_num), model="BM"))
#'  
#' require(phytools)
#' treedply(td, phylosig(phy, getVector(dat, awesomeness), "lambda", test=TRUE))
#' treedply(td, phenogram(phy, getVector(dat, SVL), ftype="off"))
#' @export
treedply.treedata <- function(tdObject, ...){
  call <- substitute(...)
  env <- new.env(parent = parent.frame(), size = 1L)
  env$phy <- tdObject$phy
  env$dat <- tdObject$data
  out <- eval(call, env)
  if(is.null(out)){
    invisible()
  } else {
    return(out)
  }
}

#' A function for returning a named vector from a data frame or matrix with row names
#' @export
getVector <- function(dat, ...){
  args <- as.character(substitute(list(...)))[-1L]
  vecs <- lapply(args,function(x) dat[[x]])
  vecs <- lapply(vecs, function(x) setNames(x, rownames(dat)))
  if(length(vecs)==1){
    vecs = vecs[[1]]
  } else {names(vecs) <- args}
  return(vecs)
}

#' @export
print.treedata <- function(tdObject, ...){
  cat("$phy \n")
  print(tdObject$phy)
  
  cat("\n$data \n")
  print(tdObject$dat)
}

