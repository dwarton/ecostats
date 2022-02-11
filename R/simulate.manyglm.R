#' Simulate from manyglm objects
#'
#' Simulates new responses for a manyglm object.
#'
#' @param formula a \code{manyglm} object from the \code{mvabund} package.
#' @param ... additional optional arguments.
#'
#' @details
#' Returns a data frame containing the response and predictors. This function just calls \code{simulate.cord} 
#' from the \code{ecoCopula} package, on a \code{cord} object constructed under default settings -- that is, 
#' it fits a copula latent variable model with two latent variables, then uses this to simulate new data. 
#' 
#' @return Simulates a data frame of new values for responses. If multiple datasets are requested, these are
#' stacked one under the other (see \code{example(simulate.cord)}. 
#' 
#' @author David Warton <david.warton@@unsw.edu.au>
#' 
#' @seealso \code{\link[mvabund]{manyglm}}, \code{\link[ecoCopula]{cord}}, \code{\link[ecoCopula]{simulate.cord}}
#'
#' @method simulate manyglm
#' @importFrom stats simulate
#' @importFrom MASS rnegbin
#' @export
simulate.manyglm = function (object, nsim=1, seed=NULL, newdata=object$data, ...) 
{
  if(ncol(fitted(object))>1)
  {
    cordObject = ecoCopula::cord(object)
    simLong    = ecoCopula:::simulate.cord(cordObject,nsim=nsim,seed=seed,newdata=newdata)
    if(nsim>1)
    {
      nDims=dim(object$residuals)
      simData = array(simLong,c(nDims[1],nsim,nDims[2]))
      simData = aperm(simData,c(1,3,2))
    }
    else
      simData=simLong
  }
  else #refit as a glm and use simulate from there
  {
    fam = switch(object$family,
                 "binomial(link=logit)"=binomial(),
                 "binomial(link=cloglog)"=binomial("cloglog"),
                 "poisson"=poisson(),
                 "gaussian"=gaussian(),
                 "gamma"=Gamma(link='log'),
                 "negative.binomial"=MASS::negative.binomial(theta=object$theta[iVar]))    
    newdata$offs = model.offset(model.frame(object$formula,data=newdata))
    if(is.null(newdata$offs))
      object.i = glm(object$formula, family=fam, data=newdata)
    else
      object.i = glm(object$formula, family=fam, data=newdata, offset=offs)
    object.i$coefficients = coef(object)
    simData = simulate(object.i,nsim=nsim,seed=seed)
  }
  
  return(simData)
}