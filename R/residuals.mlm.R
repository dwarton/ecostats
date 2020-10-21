#' Extract Multivariate Linear Model Residuals
#'
#' Residuals for a Multivariate Linear Model (\code{mlm}) object. Returns the usual (marginal) 
#' residuals by default, but optionally returns residuals from a full conditional model 
#' (\code{conditional=TRUE}), which can be used to diagnose the multivariate normality 
#' assumption.
#'
#' @param object a \code{mlm} object, typically the result of calling \code{lm} with a matrix response.
#' @param conditional logical defaults to \code{FALSE} for the usual marginal residuals, but
#' for \code{TRUE} this returns residuals from a linear model of each response against all
#' other responses, together with other model predictors.
#' @param standardize logical defaults to \code{FALSE}, but for \code{TRUE} residuals will be
#' studentized using \code{\link{rstandard}} so they are comparable across responses.
#' @param ... further arguments passed to \code{\link{residuals.lm}} or \code{\link{rstandard}}.
#'
#' @details
#' A \code{residuals} function for \code{mlm} objects, which returns the usual (marginal) 
#' residuals by default, but optionally returns residuals from a full conditional model 
#' (\code{conditional=TRUE}), that is, a linear model of each response against all responses
#' as well as predictors, which can be used to diagnose the multivariate normality assumption.
#' These can be standardized (\code{standardize=TRUE}) to facilitate overlay plots of multiple
#' responses, as is the default behaviour in \code{\link{plot.mlm}}.
#' 
#' @return A matrix of residuals
#' 
#' @author David Warton <david.warton@@unsw.edu.au>
#' 
#' @seealso \code{\link{lm}}, \code{\link{plot.mlm}}, \code{\link{predict.mlm}}, \code{\link{rstandard}}
#' @examples
#' data(iris)
#' # fit a mlm:
#' iris.mlm=lm(cbind(Sepal.Length,Sepal.Width,Petal.Length,Petal.Width)~Species,data=iris)
#' # construct full conditional residuals:
#' residuals(iris.mlm, conditional=TRUE)
#'
#' @export
residuals.mlm=function(object, conditional=FALSE, standardize = FALSE, ...)
{ 
  mf = model.frame(object)
  Ys = model.response(mf) #should be the response matrix
  X = model.matrix(object) # gets the design matrix
  nYs = dim(Ys)[2]
  resFunction = if(standardize==TRUE) rstandard else residuals
  # How do we get the response matrix from an lm object?
  # can make res = residuals if false else res = rstandard
  if(nYs==1 | conditional==FALSE) #or
    res = residuals.lm(object)
  else
  {
    res = matrix(NA,dim(Ys)[1],dim(Ys)[2])
    dimnames(res) = dimnames(Ys)
    # construct data frame with everything
    XandY = data.frame(Ys,X)
    for(iVar in 1:nYs){
      y = Ys[,iVar] #extract current response
      res[,iVar]=resFunction(lm(y~0+.,data=XandY[,-iVar]), ...)
    }
  }
  return(res)
}
