#' Add a scatterplot with a GAM smoother and confidence bands
#'
#' A scatterplot with a GAM smoother and confidence bands.
#'
#' @param x the \code{x} coordinates for the plot.
#' @param y the \code{y} coordinates for the plot.
#' @param conf.level the confidence level to use in constructing the point-wise 
#'  confidence bands.
#' @param col the colors for points. Multiple colors can be specified so that each point
#'  can be given its own color. If there are fewer colors than points they are recycled
#'  in the standard fashion.
#' @param envelope.col a character vector of length two whose first and second values 
#' respectively specify the color of the smoother and the color of the confidence band.
#' @param ... other graphical parameters passed through to \code{\link{points}}  
#'   
#'
#' @details
#' A scatterplot function which adds a GAM smoother and confidence band to the plot.
#' This is especially useful as a diagnostic tool for residual vs fits plots. The smoother
#' is fitted using \code{\link[mgcv]{gam}} and point-wise confidence bands use Wald intervals
#' (adding and subtracting a \code{qnorm(1-(1-conf.level)/2)} standard errors from predictions).
#' 
#' @return A list with the following components:
#' \item{xSorted}{A vector of \code{x} coordinates sorted from smallest to largest}
#' \item{prY    }{Predicted values for \code{y} from the GAM smoother, as a function of \code{xSorted}}
#' \item{prLow  }{Lower confidence band on the GAM smoother}
#' \item{prHi   }{Upper confidence band on the GAM smoother}
#' \item{gam.yx }{The fitted \code{gam} object}
#' 
#' @author David Warton <david.warton@@unsw.edu.au>
#' 
#' @seealso \code{\link{lm}}, \code{\link{plot.mlm}}, \code{\link{predict.mlm}}, \code{\link{rstandard}}
#' @examples
#' data(iris)
#' # fit a mlm:
#' iris.mlm=lm(cbind(Sepal.Length,Sepal.Width,Petal.Length,Petal.Width)~Species,data=iris)
#' # construct a residual vs fits plot, with GAM smoother and confidence bands:
#' plot(iris.mlm,which=1,panel=gamenvelope)
#'
#' @export
#' @importFrom mgcv gam predict.gam

gamenvelope=function(x, y, conf.level=0.95, col = par("col"), 
              envelope.col=c("blue",rgb(0,0,1,0.25)), ...)
{
  # if required, coerce x and y into vectors 
  Y = c(y)
  X = c(x)
  if(length(X)<length(Y))
  {
    X=rep(c(x),length.out=length(Y))
    warning("length of x is shorter than that of y, recycling x values")
  }    
  
  # fit a smoother using gam from mgcv package, to data sorted by X
  require(mgcv)
  xSort = sort(X,index.return=TRUE)
  gam.yx=mgcv::gam(Y[xSort$ix]~xSort$x)
  
  # make predictions and add to plot
  pr.y = mgcv::predict.gam(gam.yx,se.fit=TRUE)
  n.obs=length(xSort$ix)
  zAlpha = qnorm(1-(1-conf.level)/2)
  prHi = pr.y$fit+zAlpha*pr.y$se.fit
  prLow = pr.y$fit-zAlpha*pr.y$se.fit
  polygon(xSort$x[c(1:n.obs,n.obs:1)],c(prHi,prLow[n.obs:1]),col=envelope.col[2],border=NA) 
  lines(xSort$x,pr.y$fit,col=envelope.col[1])
  points(x,y,col=col,...)
  invisible(list(xSorted=xSort$x,prY=pr.y$fit,prLow=prLow,prHi=prHi,gam.yx=gam.yx))
}
