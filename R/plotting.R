#' Function for plotting contrast values
#' 
#' @param tree A tree of class 'phylo'
#' @param dat A named vector of continuous phenotypic data
#' @param cex.tip.label The character expansion factor for tip labels
#' @param label.offset The relative amount the tip labels are offset from the tips of the tree
#' @param ... Additional arguments passed to plot.phylo
#' 
#' @details Phylogenetic independent contrasts are calculated at each node and plotted on the phylogeny,
#' the value of of the absolute contrast is converted to a color scale. In addition, the relationship between 
#' node age and contrast value is plotted and a linear fit to the relationship between absolute contrast
#' and the age of the node (The node-height test). Under BM the relationship should be 0. 
#' @export

plotContrasts <- function(tree, dat, cex.tip.label = 0.5, label.offset=0.02, legend=TRUE, ...){
  tree$edge.length <- tree$edge.length/max(branching.times(tree))
  if(any(tree$edge.length==0)){
    tree$edge.length[tree$edge.length==0] <- min(tree$edge.length[tree$edge.length >0])/100
  }
  pic.dat <- pic(dat, tree)
  pic.col <- .color.scale(abs(pic.dat),c(0,1,1),c(1,1,0),0)
  old.par <- par()
  par(mfrow=c(2,1))
  par(mar=c(0,0,0,0))
  plot(tree, cex=cex.tip.label, label.offset=label.offset*max(branching.times(tree)), ...)
  nodelabels(pch=21, bg=pic.col)
  if (legend == TRUE) {
    pic.range <- seq(min(abs(pic.dat)), max(abs(pic.dat)), length.out=100)
    legend_image <- as.raster(matrix(rev(.color.scale(pic.range,c(0,1,1),c(1,1,0),0)), 
                                     ncol = 1))
    text(x = 1.5, y = round(seq(pic.range[1], pic.range[length(pic.range)], 
                                l = 5), 2), labels = round(seq(pic.range[1], pic.range[length(pic.range)], 
                                                               l = 5), 2))
    rasterImage(legend_image, 1.2, 0.25 * length(tree$tip.label), 
                1.21, length(tree$tip.label) - 0.25 * length(tree$tip.label))
    text(1.25, seq(0.25 * length(tree$tip.label), length(tree$tip.label) - 
                     0.25 * length(tree$tip.label), length.out = 5), 
         labels = round(seq(pic.range[1], pic.range[length(pic.range)], 
                            l = 5), 2), cex = cex.tip.label)
  }
  
  par(mar=old.par$mar)
  plot(branching.times(tree), abs(pic.dat), pch=21, xlim=c(1,0), bg=pic.col, xlab="Node Age", ylab="Absolute contrast")
  nh.lm <- lm(abs(pic.dat)~branching.times(tree))
  abline(nh.lm, lty=2, lwd=2)
}

## Utility function for plotContrasts, taken from plotrix
.color.scale <- function (x, cs1 = c(0, 1), cs2 = c(0, 1), cs3 = c(0, 1), alpha = 1, 
          extremes = NA, na.color = NA, xrange = NULL, color.spec = "rgb") {
  naxs <- is.na(x)
  if (!is.na(extremes[1])) {
    colmat <- col2rgb(extremes)
    cs1 <- colmat[1, ]/255
    cs2 <- colmat[2, ]/255
    cs3 <- colmat[3, ]/255
    color_spec <- "rgb"
  }
  maxcs1 <- ifelse(color.spec == "hcl", 360, 1)
  maxcs2 <- ifelse(color.spec == "hcl", 100, 1)
  maxcs3 <- ifelse(color.spec == "hcl", 100, 1)
  ncolors <- length(x)
  if (is.null(xrange)) {
    xrange <- range(x, na.rm = TRUE)
    drop.extremes <- FALSE
  }
  else {
    if (xrange[1] > min(x, na.rm = TRUE) || xrange[2] < max(x, 
                                                            na.rm = TRUE)) 
      stop("An explicit range for x must include the range of x values.")
    x <- c(xrange, x)
    drop.extremes = TRUE
  }
  ncs1 <- length(cs1)
  if (ncs1 > 1) {
    cs1s <- rep(cs1[ncs1], ncolors)
    xstart <- xrange[1]
    xinc <- diff(xrange)/(ncs1 - 1)
    for (seg in 1:(ncs1 - 1)) {
      segindex <- which((x >= xstart) & (x <= (xstart + 
                                                 xinc)))
      cs1s[segindex] <- .color.rescale(x[segindex], cs1[c(seg, 
                                                   seg + 1)])
      xstart <- xstart + xinc
    }
    if (min(cs1s) < 0 || max(cs1s) > maxcs1) 
      cs1s <- .color.rescale(cs1s, c(0, maxcs1))
  }
  else cs1s <- rep(cs1, ncolors)
  ncs2 <- length(cs2)
  if (ncs2 > 1) {
    cs2s <- rep(cs2[ncs2], ncolors)
    xstart <- xrange[1]
    xinc <- diff(xrange)/(ncs2 - 1)
    for (seg in 1:(ncs2 - 1)) {
      segindex <- which((x >= xstart) & (x <= (xstart + 
                                                 xinc)))
      cs2s[segindex] <- .color.rescale(x[segindex], cs2[c(seg, 
                                                   seg + 1)])
      xstart <- xstart + xinc
    }
    if (min(cs2s) < 0 || max(cs2s) > maxcs2) 
      cs2s <- .color.rescale(cs2s, c(0, maxcs2))
  }
  else cs2s <- rep(cs2, ncolors)
  ncs3 <- length(cs3)
  if (ncs3 > 1) {
    cs3s <- rep(cs3[ncs3], ncolors)
    xstart <- xrange[1]
    xinc <- diff(xrange)/(ncs3 - 1)
    for (seg in 1:(ncs3 - 1)) {
      segindex <- which((x >= xstart) & (x <= (xstart + 
                                                 xinc)))
      cs3s[segindex] <- .color.rescale(x[segindex], cs3[c(seg, 
                                                   seg + 1)])
      xstart <- xstart + xinc
    }
    if (min(cs3s) < 0 || max(cs3s) > maxcs3) 
      cs3s <- .color.rescale(cs3s, c(0, maxcs3))
  }
  else cs3s <- rep(cs3, ncolors)
  if (drop.extremes) {
    cs1s <- cs1s[-(1:2)]
    cs2s <- cs2s[-(1:2)]
    cs3s <- cs3s[-(1:2)]
  }
  xdim <- dim(x)
  colors <- do.call(color.spec, list(cs1s, cs2s, cs3s, alpha = alpha))
  if (!is.null(xdim)) 
    colors <- matrix(colors, nrow = xdim[1])
  if (length(naxs)) 
    colors[naxs] <- na.color
  return(colors)
}
## Utility function for plotContrasts, taken from plotrix
.color.rescale <- function (x, newrange) {
  if (missing(x) | missing(newrange)) {
    usage.string <- paste("Usage: rescale(x,newrange)\n", 
                          "\twhere x is a numeric object and newrange is the new min and max\n", 
                          sep = "", collapse = "")
    stop(usage.string)
  }
  if (is.numeric(x) && is.numeric(newrange)) {
    xna <- is.na(x)
    if (all(xna)) 
      return(x)
    if (any(xna)) 
      xrange <- range(x[!xna])
    else xrange <- range(x)
    if (xrange[1] == xrange[2]) 
      return(x)
    mfac <- (newrange[2] - newrange[1])/(xrange[2] - xrange[1])
    return(newrange[1] + (x - xrange[1]) * mfac)
  }
  else {
    warning("Only numeric objects can be rescaled")
    return(x)
  }
}