#' Function for checking whether a treedata object contains only numeric columns and for forcing data columns into numeric format
#'
#' This function can be used to check if a treedata object contains numeric columns and, if desired, drop all non-numeric columns.
#'
#' @param tdObject A "\code{treedata}" object
#' @param return.numeric If TRUE, then a treedata object with all numeric columns will be returned; non-numeric columns will be removed.
#' @return If return.numeric, then an object of class "\code{treedata}" with only numeric columns.
#' @examples
#' data(anolis)
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
    for(i in 1:ncol(tdObject$dat)){
      tdObject$dat[[i]] <- factor(tdObject$dat[[i]])
    }
    return(tdObject)
  } else {
    invisible()
  }
}
