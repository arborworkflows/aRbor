#' Function for making an object of class 'treedata'
#' 
#' This function generates an object of class 'treedata' that ensures that the ordering of tip labels and data remain intact. The object can be manipulated using \code{dplyr} functions.
#' 
#' @param tree An object of class 'phylo'
#' @param data A data frame or matrix
#' @param name_column The column of \code{data} that contains the names to be matched to the tree. By default, it is set to "detect" which finds the column with the most matches to the tree (including the rownames).
#' @return An object of class "\code{treedata}". The tree is pruned of tips not represented in the data, and the data is filtered for taxa not in the tree. The data is returned as a data frame tble that is compatible with \code{dplyr} functions. 
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat, name_column=1)
#' @export
make.treedata <- function(tree, data, name_column="detect") {
  if(class(tree)!="phylo") stop("tree must be of class 'phylo'")
  if(is.vector(data)){
    data <- as.matrix(data)
    colnames(data) <- "trait"
  }
  coln <- colnames(data)
  if(name_column=="detect"){
    matches <- sapply(data.frame(rownames(data), data), function(x) sum(x %in% tree$tip.label))
    if(all(matches==0)) stop("No matching names found between data and tree")
    name_column <- which(matches==max(matches))-1
  }
  if(name_column==0){
    dat.label <- rownames(data)
  } 
  dat <- tbl_df(as.data.frame(lapply(1:ncol(data), function(x) type.convert(as.character(data[,x])))))
  colnames(dat) <- coln
  #dat <- apply(dat, 2, type.convert)
  if(name_column==0){
    clnm <- colnames(dat)
    dat <- dat[,clnm, drop=FALSE]
  } else {
    if(is.numeric(name_column)){
      clnm <- (1:ncol(data))[-name_column]
    } else {
      clnm <- colnames(dat)[-which(colnames(dat)==name_column)]
    }
    dat <- dat[, clnm, drop=FALSE]
    dat.label <- as.character(data[,name_column])
  }
  data_not_tree <- setdiff(dat.label, tree$tip.label)
  tree_not_data <- setdiff(tree$tip.label, dat.label)
  phy <- drop.tip(tree, tree_not_data)
  dat <- filter(dat, dat.label %in% phy$tip.label)
  o <- match(dat.label, phy$tip.label)
  dat <- arrange(dat, o)
  td <- list(phy=phy, dat=dat)
  class(td) <- c("treedata", "list")
  attributes(td)$tip.label <- phy$tip.label
  rownames(td$dat) <- attributes(td)$tip.label
  return(td)
}

#' @export
mutate.treedata <- function(tdObject, ...){
  if(is.null(list(substitute(...))[[1]])) stop("No expressions provided to add to the treedata object")
  dat <- dplyr:::mutate_impl(tdObject$dat, dplyr:::named_dots(...), environment())
  tdObject$dat <- dat
  rownames(tdObject$dat) <- attributes(tdObject)$tip.label
  return(tdObject)
}

#' @export
select.treedata <- function(tdObject, ...){
  if(is.null(list(substitute(...))[[1]]))  stop("No criteria provided for selection")
  vars <- select_vars(names(tdObject$dat), ..., env = parent.frame())
  dat <- dplyr:::select_impl(tdObject$dat, vars)
  #dat <- lapply(dat, function(x){names(x) <- tdObject$dat[,1]; x})
  tdObject$dat <- dat
  rownames(tdObject$dat) <- attributes(tdObject)$tip.label
  return(tdObject)
}

#' @export
filter.treedata <- function(tdObject, ...){
  if(is.null(list(substitute(...))[[1]]))  stop("No criteria provided for filtering")
  tdObject$dat <- mutate(tdObject$dat, tip.label=attributes(tdObject)$tip.label)
  tdObject$dat <- filter(tdObject$dat, ...)
  attributes(tdObject)$tip.label <- tdObject$dat$tip.label
  tdObject$dat <- select(tdObject$dat, 1:(ncol(tdObject$dat)-1))
  rownames(tdObject$dat) <- attributes(tdObject)$tip.label
  tdObject$phy <- drop.tip(tdObject$phy, tdObject$phy$tip.label[!(tdObject$phy$tip.label %in% attributes(tdObject)$tip.label)])
  return(tdObject)
}

#' @export
summarize.treedata <- function(tdObject, ...){
  if(is.null(list(substitute(...))[[1]]))  stop("No expression provided to summarize data")
  res <- summarize(tdObject$dat, ...)
  return(res)
}

#' @export
summarise.treedata <- function(tdObject, ...){
  if(is.null(list(substitute(...))[[1]])) stop("No expression provided to summarize data")
  res <- summarise(tdObject$dat, ...)
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

#' @rdname treeply.treedata
#' @export
treeply <- function(tdObject, ...){
  UseMethod("treeply")
}

#' Run a function on the phylogeny of a 'treedata' object
#' @description Applies a function to the phylogeny in a 'treedata' object. If the order of tips are changed, or if tips are dropped, then the data are automatically reordered to match the tree.
#' @param tdObject An object of class 'treedata'
#' @param FUN A function that operates on an object of class 'phylo'
#' 
#' @return An object of class 'treedata'
#' 
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat, name_column=1)
#' td_OU <- treeply(td, rescale, model="OU", 10)
#'   
#' par(mfrow=c(1,2))
#' plot(td$phy)
#' plot(td_OU$phy)
#' 
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

#' @rdname treedply.treedata
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
  if(!is.call(substitute(...))){
    call <- list(...)[[1]]
  } else {
    call <- substitute(...)
  }
  env <- new.env(parent = parent.frame(), size = 1L)
  env$phy <- tdObject$phy
  env$dat <- tdObject$dat
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
  arg_sub <- type.convert(args)
  if(is.numeric(arg_sub) | is.integer(arg_sub)) args <- arg_sub
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

#' @export
checkNumeric <- function(tdObject, return.numeric=TRUE) {
  valid <- which(sapply(tdObject$dat, is.numeric))
  if(length(valid) < ncol(tdObject$dat)){
    if(length(valid)==0){
      stop("Dataset does not contain any numeric data that can be used for continuous ancestral state reconstruction") }
    else {
      not.valid <- colnames(tdObject$dat)[which(sapply(tdObject$dat, is.numeric))]
      warning(paste("Not all data continuous, dropping non-numeric data columns:", paste(not.valid, collapse=" ")))
      tdObject <- select(tdObject, valid)
    }
  }
  if(return.numeric){
    return(tdObject)
  } else {
    invisible()
  }
}

#' @export
checkFactor <- function(tdObject, return.factor=TRUE) {
  classes <- sapply(tdObject$dat, class)
  valid <- which(classes=="factor")
  if(length(valid) < ncol(tdObject$dat)){
    #Which data are numeric
    are.numeric <- which(classes %in% c("numeric", "integer"))
    if(length(are.numeric) > 0){
      warning("Data contain numeric entries, which will be converted to factors for discrete ancestral state reconstruction")
      #convert them to factors
      tdObject$dat[, are.numeric] <- lapply(tdObject$dat[,are.numeric], factor)
      ##Check to see if converted data has any columns that appear continuous
      classes <- sapply(tdObject$dat, class)
      valid <- which(classes=="factor")
      too.many.levels <- which(sapply(tdObject$dat,function(x) length(unique(x)))==nrow(tdObject$dat))
      if(length(too.many.levels) > 0){
        warning(paste("Conversion failed for data columns", paste(colnames(tdObject$dat)[too.many.levels], collapse=" "), "as these data have no shared states. These data will be removed", sep=" "))
        valid <- valid[which(!(valid %in% too.many.levels))]
      }
    }

    if(length(valid)==0){
      stop("Data does not contain any data that can be used for discrete ancestral state reconstruction") 
    } else {
      tdObject$dat <- select(tdObject$dat, valid)
    }
  }
  if(return.factor){
    return(tdObject)
  } else {
    invisible()
  }
}

