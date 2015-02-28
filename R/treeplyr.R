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
  if(is.null(colnames(data))){
    colnames(data) <- paste("trait", 1:ncol(data), sep="")
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
  dat <- tbl_df(as.data.frame(lapply(1:ncol(data), function(x) type.convert(apply(data[,x, drop=FALSE], 1, as.character)))))
  #dat <- tbl_df(as.data.frame(lapply(1:ncol(data), function(x) type.convert(as.character(data[,x])))))
  #dat <- data
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
    dat.label <- as.character(data[[name_column]])
  }
  data_not_tree <- setdiff(dat.label, tree$tip.label)
  tree_not_data <- setdiff(tree$tip.label, dat.label)
  phy <- drop.tip(tree, tree_not_data)
  dat <- filter(dat, dat.label %in% phy$tip.label)
  dat.label <- dat.label[dat.label %in% phy$tip.label]
  if(any(duplicated(dat.label))){
    warning("Duplicated data in dataset, selecting first unique entry for each species")
    dat <- filter(dat, !duplicated(dat.label))
    dat.label <- dat.label[!duplicated(dat.label)]
  }
  o <- match(dat.label, phy$tip.label)
  dat <- arrange(dat, o)
  td <- list(phy=phy, dat=dat)
  class(td) <- c("treedata", "list")
  attributes(td)$tip.label <- phy$tip.label
  rownames(td$dat) <- attributes(td)$tip.label
  return(td)
}

#' Function for mutating an object of class 'treedata'
#' 
#' This function can be used to add new variables to a treedata object; see \code{\link{mutate}}.
#' 
#' @param tdObject A "\code{treedata}" object
#' @return An object of class "\code{treedata}" with new data added. 
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat[,-(3:5)], name_column=1)
#' tdmutate <- mutate(td, anolis$dat[,3])
#' @export
mutate.treedata <- function(tdObject, ...){
  if(is.null(list(substitute(...))[[1]])) stop("No expressions provided to add to the treedata object")
  dat <- mutate_impl.dplyr(tdObject$dat, named_dots.dplyr(...), environment())
  tdObject$dat <- dat
  rownames(tdObject$dat) <- attributes(tdObject)$tip.label
  return(tdObject)
}

#' Function for selecting columns from an object of class 'treedata'
#' 
#' This function can be used to select a subset of variables (columns) from a treedata object; see \code{\link{select}}.
#' 
#' @param tdObject A "\code{treedata}" object
#' @return An object of class "\code{treedata}" with specified variables selected. 
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat, name_column=1)
#' tdselect <- select(td, SVL)
#' @export
select.treedata <- function(tdObject, ...){
  if(is.null(list(substitute(...))[[1]]))  stop("No criteria provided for selection")
  vars <- select_vars(names(tdObject$dat), ..., env = parent.frame())
  dat <- select_impl.dplyr(tdObject$dat, vars)
  #dat <- lapply(dat, function(x){names(x) <- tdObject$dat[,1]; x})
  tdObject$dat <- dat
  rownames(tdObject$dat) <- attributes(tdObject)$tip.label
  return(tdObject)
}

#' Function for filtering rows from an object of class 'treedata'
#' 
#' This function can be used to select a subset of species (rows) from a treedata object; see \code{\link{filter}}.
#' 
#' @param tdObject A "\code{treedata}" object
#' @return An object of class "\code{treedata}" with specified species selected. 
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat, name_column=1)
#' tdfilter <- filter(td, island=="Cuba")
#' @export
filter.treedata <- function(tdObject, ...){
  #if(is.null(list(substitute(...))[[1]]))  stop("No criteria provided for filtering")
  tdObject$dat <- mutate(tdObject$dat, tip.label=attributes(tdObject)$tip.label)
  tdObject$dat <- filter(tdObject$dat, ...)
  attributes(tdObject)$tip.label <- tdObject$dat$tip.label
  tdObject$dat <- select(tdObject$dat, 1:(ncol(tdObject$dat)-1))
  rownames(tdObject$dat) <- attributes(tdObject)$tip.label
  tdObject$phy <- drop.tip(tdObject$phy, tdObject$phy$tip.label[!(tdObject$phy$tip.label %in% attributes(tdObject)$tip.label)])
  return(tdObject)
}

#' @name summarise
#' @aliases summarize summarize.treedata summarise.treedata
#' @title Function for summarizing an object of class 'treedata'
#' 
#' This function can be used to summarize a treedata object.
#' 
#' @param tdObject A "\code{treedata}" object
#' @return A summary of the treedata object 
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat, name_column=1)
#' summarize(td)
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

#' @rdname reorder.treedata
#' @export
reorder <- function(tdObject, ...){
  UseMethod("reorder")
}

#' Reorder a 'treedata' object
#' @description Reorders a 'treedata' object. Both the tips and the data are automatically reordered to match.
#' @param tdObject An object of class 'treedata'
#' @param order Method for reordering
#' @param index.only XXX
#' 
#' @return An object of class 'treedata'
#' 
#' @examples
#' x <- 10
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
    tdObject$dat <- tdObject$dat[match(tdObject$phy$tip.label, attributes(tdObject)$tip.label),]
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
#' treedply(td, fitContinuous(phy, getVector(dat, SVL), model="BM", ncores=2))
#' #treedply(td, fitContinuous(phy, select(dat, PCI_limbs, PCII_head, 
#'   #PCIII_padwidth_vs_tail, PCIV_lamella_num), model="BM", ncores=2))
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
print.treedata <- function(x, ...){
  cat("$phy \n")
  print(x$phy)
  
  cat("\n$data \n")
  print(x$dat)
}


#' Function for checking whether a treedata object contains only numeric columns and for forcing data columns into numeric format
#' 
#' This function can be used to check if a treedata object contains numeric columns and, if desired, drop all non-numeric columns.
#' 
#' @param tdObject A "\code{treedata}" object
#' @param return.numeric If TRUE, then a treedata object with all numeric columns will be returned; non-numeric columns will be removed.
#' @return If return.numeric, then an object of class "\code{treedata}" with only numeric columns. 
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat, name_column=1)
#' tdnumeric <- checkNumeric(td)
#' @export
checkNumeric <- function(tdObject, return.numeric=TRUE) {
  valid <- which(sapply(tdObject$dat, is.numeric))
  if(length(valid) < ncol(tdObject$dat)){
    if(length(valid)==0){
      stop("Dataset does not contain any numeric data that can be used for continuous ancestral state reconstruction") }
    else {
      not.valid <- colnames(tdObject$dat)[which(!sapply(tdObject$dat, is.numeric))]
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

#' Function for checking whether a treedata object contains only factors and for forcing data columns into factor format
#' 
#' This function can be used to check if a treedata object contains factors and, if desired, convert all columns automatically to factors.
#' 
#' @param tdObject A "\code{treedata}" object
#' @param return.factor If TRUE, then a treedata object with all factors will be returned; columns will be forced into factors using \code{factor} and any with no repeated elements will be removed.
#' @return If return.factor, then an object of class "\code{treedata}" with all columns as factors. 
#' @examples
#' data(anolis)
#' td <- make.treedata(anolis$phy, anolis$dat, name_column=1)
#' tdforcefactor <- checkFactor(td)
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
      tdObject$dat[, are.numeric] <- lapply(as.data.frame(tdObject$dat[,are.numeric]), factor)
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

# Helper functions that allow string arguments for  dplyr's data modification functions like arrange, select etc. 
# Author: Sebastian Kranz

# Examples are below

#' Modified version of dplyr's filter that uses string arguments
#' @export
s_filter = function(.data, ...) {
  eval.string.dplyr(.data,"filter", ...)
}

#' Modified version of dplyr's select that uses string arguments
#' @export
s_select = function(.data, ...) {
  eval.string.dplyr(.data,"select", ...)
}

#' Modified version of dplyr's arrange that uses string arguments
#' @export
s_arrange = function(.data, ...) {
  eval.string.dplyr(.data,"arrange", ...)
}

#' Modified version of dplyr's arrange that uses string arguments
#' @export
s_mutate = function(.data, ...) {
  eval.string.dplyr(.data,"mutate", ...)
}

#' Modified version of dplyr's summarise that uses string arguments
#' @export
s_summarise = function(.data, ...) {
  eval.string.dplyr(.data,"summarise", ...)
}

#' Modified version of dplyr's group_by that uses string arguments
#' @export
s_group_by = function(.data, ...) {
  eval.string.dplyr(.data,"group_by", ...)
}

#' Modified version of dplyr's group_by that uses string arguments
#' @export
s_treeply = function(.data, ...) {
  eval.string.dplyr(.data,"treeply", ...)
}

#' Modified version of dplyr's group_by that uses string arguments
#' @export
s_treedply = function(.data, ...) {
  eval.string.dplyr(.data,"treedply", ...)
}

#' Internal function used by s_filter, s_select etc.
eval.string.dplyr = function(.data, .fun.name, ...) {
  args = list(...)
  args = unlist(args)
  code = paste0(.fun.name,"(.data,", paste0(args, collapse=","), ")")
  df = eval(parse(text=code,srcfile=NULL))
  df  
}

#' @export
select_.treedata <- function(.data, ..., .dots){
  dots <- all_dots(.dots, ...)
  vars <- select_vars_(names(.data$dat), dots)
  dat <- .data$dat[, vars, drop = FALSE]
  row.names(dat) <- attributes(.data)$tip.label
  .data$dat <- dat
  return(.data)
}

#' @export
mutate_.treedata <- function(.data, ..., .dots){
  dots <- all_dots(.dots, ..., all_named=TRUE)
  dat <- mutate_impl.dplyr(.data$dat, dots)
  row.names(dat) <- attributes(.data)$tip.label
  .data$dat <- dat
  return(.data)
}

#' @export
filter_.treedata <- function(.data, ..., .dots){
  dots <- all_dots(.dots, ..., all_named=TRUE)
  .data$dat <- mutate(.data$dat, tip.label=attributes(.data)$tip.label)
  dat <- filter_impl.dplyr(.data$dat, dots)
  .data$dat <- dat
  attributes(.data)$tip.label <- .data$dat$tip.label
  .data$dat <- select(.data$dat, 1:(ncol(.data$dat)-1))
  .data$phy <- drop.tip(.data$phy, .data$phy$tip.label[!(.data$phy$tip.label %in% attributes(.data)$tip.label)])
  return(.data)
}

#' @export
select_impl.dplyr <- function (df, vars) 
{
    .Call("dplyr_select_impl", PACKAGE = "dplyr", df, vars)
}

#' @export
filter_impl.dplyr <- function (df, dots) 
{
    .Call("dplyr_filter_impl", PACKAGE = "dplyr", df, dots)
}

#' @export
mutate_impl.dplyr <- function (df, dots) 
{
    .Call("dplyr_mutate_impl", PACKAGE = "dplyr", df, dots)
}

#' @export
named_dots.dplyr <- function (...) 
{
    auto_name(dots(...))
}


