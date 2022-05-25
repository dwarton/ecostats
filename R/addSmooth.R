#' addSmooth Add a smoother to a plot, with pointwise confidence band
#'
#' Adds a smoother to a plot of \code{y} against \code{x}, and a confidence band around the smoother with pointwise coverage
#' (*not* global, unlike \code{\link{plotenvelope}}). The smoother is constructed using \code{\link[mgcv]{gam}}
#' and confidence bands use Wald intervals, extracting the standard error from \code{\link[mgcv]{predict.gam}}.
#' 
#' @param x is the \code{x} co-ordinates for the plot, _e.g._ the fitted values in a plot of residuals vs fitted values.
#' @param y is the \code{y} co-ordinates for the plot
#' @param conf.level the confidence level to use in constructing the confidence band around the smoother.
#' @param line.col color of smoother. Defaults to "slateblue4".
#' @param envelope.col color of the global envelope around the expected trend. All data points should always stay within this envelope
#' (and will for a proportion \code{conf.level} of datasets satisfying model assumptions). If a smoother has been used on
#' the residual vs fits or scale-location plot, the envelope is constructed around the smoother, that is, the smoother should always stay
#' within the envelope.
#' @param ... further arguments sent through to \code{lines} and \code{polygon}.
#' 
#' @details
#' A challenge when interpreting diagnostic plots is understanding the extent to which
#' deviations from the expected pattern could be due to random noise (sampling variation)
#' rather than actual assumption violations. This function is intended to assess this, 
#' by simulating multiple realizations of residuals (and fitted values) in situations where 
#' assumptions are satisfied, and plotting a global simulation envelope around these at level \code{conf.level}.
#' 
#' This function adds a smoother to a plot of \code{y} against \code{x} using \code{\link[mgcv]{gam}}, and a confidence band 
#' around the smoother with pointwise coverage (*not* global, unlike \code{\link{plotenvelope}}) using Wald intervals
#' that take the standard error of predicted values from \code{\link[mgcv]{predict.gam}}.
#' 
#' When constructing a plot to diagnose model fit, \code{\link{plotenvelope}} is preferred, because this constructs
#' simulation envelopes globally, simplifying interpretation, and tends to give more accurate envelopes, via simulation.
#' This function is intended to fit trends to data for descriptive purposes, but it could be used (with caution) 
#' to diagnose model fit in residual plots in situations where \code{\link{plotenvelope}} is unavailable.
#' 
#' @return A smoother is added to the existing plot, with confidence band. In addition, a data frame is invisibly
#' returned with the following components:\describe{
#' \item{\code{X}}{The values of \code{x} at which the smoother is evaluated.}
#' \item{\code{Yhat}}{Values of the smoother used to predict \code{y}, for each value of \code{x}.}
#' \item{\code{Yhi}}{The upper bound on the confidence band, for each value of \code{x}.}
#' \item{\code{Ylo}}{The lower bound on the global envelope, for each value of \code{x}.}
#' }
#' 
#' @author David Warton <david.warton@@unsw.edu.au>
#' @seealso \code{\link{plotenvelope}}
#'
#' @references Warton DI (2022) Eco-Stats - Data Analysis in Ecology, from \emph{t}-tests to multivariate abundances. Springer, ISBN 978-3-030-88442-0
#'
#' @examples
#' data(globalPlants)
#' plot(log(height)~lat,data=globalPlants)
#' with(globalPlants, addSmooth(lat,log(height)) )
#' 
#' @import graphics

#' @export
addSmooth=function(x,y,conf.level=0.05,line.col="slateblue4", envelope.col = adjustcolor(line.col, alpha.f = 0.1), ...)
{
  nPoints  = 500
  dat      = data.frame(X=as.vector(x),Y=as.vector(y))
  ftGAM    = mgcv::gam(Y~X,data=dat)
  Xpred    = seq(min(dat$X),max(dat$X),length=nPoints)
  smoothDF = data.frame(X=Xpred)
  prGAM    = mgcv::predict.gam(ftGAM,newdata=smoothDF,se.fit=TRUE)

  zAlpha   = qnorm(1-conf.level/2)
  smoothDF$Yhat = prGAM$fit
  smoothDF$Yhi  = prGAM$fit + zAlpha*prGAM$se.fit
  smoothDF$Ylo  = prGAM$fit - zAlpha*prGAM$se.fit
  
  lines(Xpred,prGAM$fit,col=line.col,...)
  
  polygon( c(smoothDF$X,smoothDF$X[nPoints:1]), c(smoothDF$Yhi,smoothDF$Ylo[nPoints:1]), col=envelope.col, border=NA, ...)
  
  invisible(smoothDF)
}
