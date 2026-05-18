plotGradientSimple <- function(hM, Gradient, predY, measure, xlabel = NULL, ylabel = NULL, 
                         index = 1, q = c(0.025, 0.5, 0.975), cicol = rgb(0, 0, 1, alpha = 0.5), 
                         pointcol = "lightgrey", pointsize = 1, showData = FALSE, 
                         jigger = 0, yshow = NA, showPosteriorSupport = TRUE, main, 
                         x_ticks = NULL, x_labels = NULL, ...) {
  Pr <- NA
  if (is.null(xlabel)) {
    xlabel <- if (is.matrix(hM$X)) colnames(Gradient$XDataNew)[1] else colnames(Gradient$XDataNew[[1]])[1]
  }
  xx <- if (is.matrix(hM$X)) Gradient$XDataNew[, 1] else if (measure == "Y") Gradient$XDataNew[[index]][, 1] else Gradient$XDataNew[[1]][, 1]
  ngrid <- length(xx)
  
  if (measure == "S") {
    predS <- abind(lapply(predY, rowSums), along = 2)
    Pr <- mean(predS[ngrid, ] > predS[1, ])
    qpred <- apply(predS, 1, quantile, probs = q, na.rm = TRUE)
    ylabel <- if (is.null(ylabel)) if (all(hM$distr[, 1] == 2)) "Species richness" else if (all(hM$distr[, 1] == 3)) "Total count" else "Summed response"
  } else if (measure == "Y") {
    tmp <- abind(predY, along = 3)
    Pr <- mean(tmp[ngrid, index, ] > tmp[1, index, ])
    qpred <- apply(tmp, c(1, 2), quantile, probs = q, na.rm = TRUE)[, , index]
    ylabel <- if (is.null(ylabel)) hM$spNames[[index]] else ylabel
  } else if (measure == "T") {
    predT <- if (all(hM$distr[, 1] == 1)) lapply(predY, function(a) (exp(a) %*% hM$Tr) / matrix(rep(rowSums(exp(a)), hM$nt), ncol = hM$nt)) else lapply(predY, function(a) (a %*% hM$Tr) / matrix(rep(rowSums(a), hM$nt), ncol = hM$nt))
    predT <- abind(predT, along = 3)
    Pr <- mean(predT[ngrid, index, ] > predT[1, index, ])
    qpred <- apply(predT, c(1, 2), quantile, probs = q, na.rm = TRUE)[, , index]
    ylabel <- if (is.null(ylabel)) hM$trNames[[index]] else ylabel
  }
  
  lo <- qpred[1, ]
  hi <- qpred[3, ]
  me <- qpred[2, ]
  lo1 <- min(lo, yshow, na.rm = TRUE)
  hi1 <- max(hi, yshow, na.rm = TRUE)
  
  if (showData) {
    XDatacol <- if (is.matrix(hM$X)) which(colnames(Gradient$XDataNew)[1] == colnames(hM$XData)) else which(colnames(Gradient$XDataNew[[1]])[1] == colnames(hM$XData[[1]]))
    pY <- if (measure == "S") rowSums(hM$Y, na.rm = TRUE) else if (measure == "Y") hM$Y[, index] else if (all(hM$distr[, 1] == 1)) (exp(hM$Y) %*% hM$Tr) / matrix(rep(rowSums(exp(hM$Y)), hM$nt), ncol = hM$nt)[, index] else (hM$Y %*% hM$Tr) / matrix(rep(rowSums(hM$Y), hM$nt), ncol = hM$nt)[, index]
    pX <- if (is.matrix(hM$X)) hM$XData[, XDatacol] else hM$XData[[1]][, XDatacol]
    hi1 <- max(hi1, max(pY, na.rm = TRUE))
    lo1 <- min(lo1, min(pY, na.rm = TRUE))
  }
  
  if (!is.null(x_ticks) && !is.null(x_labels)) {
    plot(xx, me, ylim = c(lo1, hi1), type = "l", xlab = xlabel, ylab = ylabel, xaxt = 'n', ...)
    axis(1, at = x_ticks, labels = x_labels)
  } else {
    plot(xx, me, ylim = c(lo1, hi1), type = "l", xlab = xlabel, ylab = ylabel, ...)
  }
  
  polygon(c(xx, rev(xx)), c(lo, rev(hi)), col = cicol, border = FALSE)
  lines(xx, me, lwd = 2)
  
  if (showData) {
    if (jigger > 0) {
      de <- (hi1 - lo1) * jigger
      pY <- lo1 + de + (hi1 - lo1 - 2 * de) * (pY - lo1) / (hi1 - lo1) + runif(length(pY), min = -jigger, max = jigger)
    }
    points(pX, pY, pch = 16, col = pointcol, cex = pointsize)
  }
  
  if (!missing(main)) title(main = main)
  if (showPosteriorSupport) mtext(gettextf("Pr[pred(%s=min) %s pred(%s=max)] = %.2f", xlabel, ifelse(Pr < 0.5, ">", "<"), xlabel, ifelse(Pr < 0.5, 1 - Pr, Pr)))
  
  Pr
}
