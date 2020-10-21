#' Predict method for Multivariate Linear Model Fits
#'
#' Predicted values based on a multivariate linear model (\code{mlm}) object.
#' Returns the usual (marginal) predictions by default, but optionally returns
#' predictions from a full conditional model (\code{conditional=TRUE}), which
#' can be used to diagnose the multivariate normality assumption.
#'
#' @param object a \code{mlm} object, typically the result of calling \code{lm} with a matrix response.
#' @param conditional logical defaults to \code{FALSE} for the usual marginal predictions, but
#' for \code{TRUE} this returns predictions for each respones conditional on the values of
#' all other responses.
#' @param standardize logical defaults to \code{FALSE}, but for \code{TRUE} predictions will be
#' standardised so they are comparable across responses.
#' @param ... further arguments passed to \code{predict.lm}.
#'
#' @details
#' A \code{predict} function for \code{mlm} objects, which returns the usual (marginal) 
#' predicted values by default, but optionally returns predictions from a full conditional 
#' model (\code{conditional=TRUE}), that is, from a linear model for each response as a function
#' of all responses as well as predictors. This can be used in plots to diagnose the
#' multivariate normality assumption.
#' These can be standardized (\code{standardize=TRUE}) to facilitate overlay plots of multiple
#' responses, as is the default behaviour in \code{\link{plot.mlm}}.
#' 
#' @return A matrix of predicted values
#' 
#' @author David Warton <david.warton@@unsw.edu.au>
#' 
#' @seealso \code{\link{lm}}, \code{\link{plot.mlm}}, \code{\link{residuals.mlm}}
#' @examples
#' data(iris)
#' # fit a mlm:
#' iris.mlm=lm(cbind(Sepal.Length,Sepal.Width,Petal.Length,Petal.Width)~Species,data=iris)
#' # predict each response conditionally on the values of all other responses:
#' predict(iris.mlm, conditional=TRUE)
#' 
#' @method predict mlm
#' @S3method predict mlm
#' @export predict.mlm
predict.mlm = function(object, conditional = FALSE, standardize=FALSE, ...)
{
  mf = model.frame(object) 
  Ys = model.response(mf)
  Ys = as.matrix(Ys) # think this fixes dim = NULL problemsfit
  X = model.matrix(object) # gets the design matrix
  nYs = dim(Ys)[2]
  
  results = matrix(NA,dim(Ys)[1],nYs)
  dimnames(results) = dimnames(Ys)
  if(nYs == 1 | conditional == FALSE)
    for (iVar in 1:nYs){
      y = Ys[ ,iVar]
      fit = lm(y~X)
      results[ ,iVar] = predict(fit, ...)
    }
  else
  {
    XandY = data.frame(Ys, X)
    for (iVar in 1:nYs){
      y = Ys[ ,iVar]
      fit = lm(y~.,data=XandY[,-iVar])
      results[ ,iVar] = predict(fit, ...)
    }
  }
  if(standardize==TRUE)
    results=scale(results)
  return(results)
}
