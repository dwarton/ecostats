#' Simulate Responses from a Multivariate Linear Model
#'
#' Simulate one or more sets of responses from  a Multivariate Linear Model (\code{mlm}) object.
#'
#' @param object a \code{mlm} object, typically the result of calling \code{lm} with a matrix response.
#' @param nsim number of replicate datasets to simulate. If \code{nsim} is greater than 1
#' the output is arranged in a 3D array.
#' @param seed an object specifying if and how the random number generator should be 
#' initialized (‘seeded’). Either NULL or an integer that will be used in a call to set.seed
#' before simulating the response vectors. If set, the value is saved as the "seed" attribute 
#' of the returned value. The default, NULL will not change the random generator state, and 
#' return .Random.seed as the "seed" attribute, see ‘Value’
#' @param ... additional optional arguments.
#'
#' @details
#' A \code{simulate} function for \code{mlm} objects, which simulates one or more sets of responses from 
#' a Multivariate Linear Model (\code{mlm}) object. If multiple sets of responses are requested,
#' they are returned in a 3D array, with simulation number along the third dimension.
#' Weights argument is currently ignored -- constant variance-covariance matrix assumed.
#' 
#' @return A matrix of simulated values for the response
#' 
#' @author David Warton <david.warton@@unsw.edu.au>
#' 
#' @seealso \code{\link{lm}}, \code{\link{plot.mlm}}, \code{\link{predict.mlm}}, \code{\link{rstandard}}
#' @examples
#' data(iris)
#' # fit a mlm:
#' iris.mlm=lm(cbind(Sepal.Length,Sepal.Width,Petal.Length,Petal.Width)~Species,data=iris)
#' # simulate new responses:
#' simulate.mlm(iris.mlm)
#'
#' @method simulate mlm
#' @S3method simulate mlm
#' @export simulate.mlm
#' @importFrom mvtnorm rmvnorm

#simulate function needed for mlm objects
simulate.mlm = function(object, nsim=1, seed=NULL, ...)
{
  # code chunks taken from stats:::simulate.lm R version 4.0.3ish on 30/7/20
  # get random seed
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
    runif(1)
  if (is.null(seed)) 
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  
  # set dimensions of simulated data object
  ftd <- fitted(object)
  nm <- dimnames(ftd)
  n <- nrow(ftd)
  vars <- var(residuals(object))
  #  if (!(is.null(w <- object$weights) || (length(w) == 
  #                                         1L && w == 1))) vars <- vars/w
  val <- replicate(nsim,ftd+mvtnorm::rmvnorm(n,sigma=vars))
  colnames(val)=nm[[2]]
  rownames(val) <- nm[[1]]
  dimnames(val)[[3]] <- paste0("sim_",1:nsim)
  attr(val, "seed") <- RNGstate
  val
}

#setMethod("simulate", "mlm", simulate.mlm)
