#' Function that uses  "leave 1 out" validation of a discrete character model
#' 
#' @param asr An asrArbor object produced by aceArbor
#' @param prior Either "stationary" or "equal" 
#' @param plot A logical indicating whether a plot should be produced
#' @param pal A color palette to color the nodes
validateASR <- function(asr, prior="stationary", plot=TRUE, pal="rainbow"){
  ntips <- length(attributes(asr)$td$phy$tip.label)
  nnodes <- attributes(asr)$td$phy$Nnode
  model <- attributes(asr)$discreteModelType
  tree <- attributes(asr)$td$phy
  states <- factor(setNames(attributes(asr)$td$dat[,1], tree$tip.label))
  fit <- attributes(asr[[1]])$fit
  names(fit)
  npars <- length(fit$par.full)
  k <- (1+sqrt(1+4*npars))/2
  Q <- matrix(0, ncol=k, nrow=k)
  Q <- sapply(1:k, function(x) {Q[x, -x] <- fit$par.full[(x-1)*(k-1)+(1:(k-1))]; Q[x, ]})
  diag(Q) <- apply(Q, 1, sum)*-1
  XX <- asr[[1]]
  attributes(XX)$fit <- NULL
  X <- matrix(0, nrow=k, ncol=ntips)
  X[(0:(ntips-1))*k + as.numeric(states)] <- 1
  if(prior=="equal"){
    pis <- rep(1/k, k)
  }
  if(prior=="stationary"){
    pis <- stationary.freqs(Q)
  }
  X <- t(X)
  rownames(X) <- names(states)
  if(prior=="equal"){
    pis <- rep(1/k, k)
  }
  if(prior=="stationary"){
    pis <- stationary.freqs(Q)
  }
  probs <- sapply(1:ntips,function(x){X[x,] <- pis; .reRoot1(x, tree, X, Q)})
  #probs <- t(probs)
  #dd <- barplot(probs, col = cols, horiz=TRUE, xlim=c(-0.4, 1))
  #points(rep(-0.05, ntips),dd, pch=22, bg=cols[states], col="transparent")
  #points(rep(-0.04, ntips),dd, pch=22, bg=cols[states], col="transparent")
  #points(rep(-0.03, ntips),dd, pch=22, bg=cols[states], col="transparent")
  #points(rep(-0.02, ntips),dd, pch=22, bg=cols[states], col="transparent")
  #text(-0.05, y=dd, labels=tree$tip.label, pos=2, cex=0.5)
  #add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
  #                  y=max(dd), fsize=0.8)
  if(plot){
    cols<-setNames(pal(k)[1:k],sort(unique(states)))
    Probs <- t(probs)
    plot(tree, type="fan", cex=cex.lab, label.offset=label.offset)#, ...)
    nodelabels(node=(ntips+1):(ntips+nnodes), pie=XX,piecol=cols,cex=cex.node)
    #tiplabels(tip=probs, pch=21, cex=2, bg="white", col="yellow")
    
    tiplabels(pch=21, bg=cols[states], col=cols[states], cex=cex.tip*5)
    tiplabels(pie=Probs,piecol=cols,cex=cex.tip)
    add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                      y=-max(nodeHeights(tree)))
  }
  
  best <- apply(t(probs), 1, function(x) which(x==max(x)))
  piCorrect <- sum(best==as.numeric(states), na.rm=TRUE)/length(states[!is.na(states)])
  correct <- tapply(best==as.numeric(states), states, sum)
  attempted <- tapply(best==as.numeric(states), states, length)
  pct <- correct/attempted
  cumulative <- sapply(1:k, function(y) sum(sapply((1:ntips)[states==levels(states)[y]], function(x) probs[as.numeric(states)[x], x])))
  exp.cs <- as.vector(pis)*tapply(states, states, length)
  prob.ratio <- cumulative/exp.cs
  res <- data.frame(correct, attempted, "pct"=round(pct,4), "prior"=round(as.vector(pis),3), "cumulativeProb"=round(cumulative,3), "expProb"=round(exp.cs,3), "ProbRatio"=round(prob.ratio,3))
  return(list(probs=Probs, summary=res, Q=Q))
}

#' A function to predict the states of missing tips
#' 
#' @param asr An asrArbor object produced by aceArbor
#' @param prior Either "stationary" or "equal" 
#' @param plot A logical indicating whether a plot should be produced
#' @param cex.node Character expansion vector for the node pie charts
#' @param cex.tips Character expansion vector for the tip pie charts
#' @param cex.lab Character expansion vector for the tip labels
#' @param pal A color palette to color the nodes
predictMissingTips <- function(asr, prior="stationary", plot=TRUE, cex.node=0.5, cex.tip=0.5, cex.lab=1, pal=rainbow, label.offset=0.1, ...){
  if(attributes(asr)$charType != "discrete") stop("Must be an asrArbor object with charType=='discrete'")
  fit <- attributes(asr[[1]])$fit
  names(fit)
  npars <- length(fit$par.full)
  k <- (1+sqrt(1+4*npars))/2
  Q <- matrix(0, ncol=k, nrow=k)
  Q <- sapply(1:k, function(x) {Q[x, -x] <- fit$par.full[(x-1)*(k-1)+(1:(k-1))]; Q[x, ]})
  diag(Q) <- apply(Q, 1, sum)*-1
  tree <- attributes(asr)$td$phy
  states <- setNames(attributes(asr)$td$dat[,1], tree$tip.label)
  missing <- which(is.na(states))
  ntips <- length(tree$tip.label)
  nnodes <- tree$Nnode
  X <- matrix(0, nrow=k, ncol=ntips)
  X[(0:(ntips-1))*k + as.numeric(states)] <- 1
  if(prior=="equal"){
    pis <- rep(1/k, k)
  }
  if(prior=="stationary"){
    pis <- stationary.freqs(Q)
  }
  X[,missing] <- pis 
  X <- t(X)
  rownames(X) <- names(states)
  model <- attributes(asr)$discreteModelType
  res <- rerootingMethod(tree, X, model, fixedQ=Q)
  cols<-setNames(pal(k)[1:k],sort(unique(states)))
  if(plot){
    plot(tree, type="fan", cex=cex.lab, label.offset=label.offset, ...)
    nodelabels(node=as.numeric(rownames(res$marginal.anc[(ntips+1):(ntips+nnodes),])),
               pie=res$marginal.anc[(ntips+1):(ntips+nnodes),],piecol=cols,cex=cex.node)
    tiplabels(tip=missing, pch=21, cex=2, bg="yellow")
    tiplabels(pie=res$marginal.anc[1:ntips,],piecol=cols,cex=cex.tip)
    add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                      y=-max(nodeHeights(tree)))
  }
  res$missing <- res$marginal.anc[missing,]
  res$Q <- Q
  invisible(res)
}

#' Internal function for calculating the marginal ancestral state of a node using the rerooting method. Adapted from
#' phytools. 
.reRoot1 <- function (i, tree, x, Q, model) {
  n <- length(tree$tip.label)
  yy <- x
  yy <- yy[tree$tip.label, ]
  yy <- yy/rowSums(yy)
  tt <- reroot(tree, node.number = i, position = tree$edge.length[which(tree$edge[,2] == i)])
  res.i <- suppressWarnings(apeAce(tt, yy, model = model, fixedQ = Q)$lik.anc[1,])
  return(res.i)
}
