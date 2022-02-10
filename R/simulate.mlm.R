#' Simulate Responses from a Multivariate Linear Model
#'
#' Simulate one or more sets of responses from a Multivariate Linear Model (\code{mlm}) object.
#'
#' @param object a \code{mlm} object, typically the result of calling \code{lm} 
#' where the response is a matrix.
#' @param nsim number of replicate datasets to simulate. If \code{nsim} is greater than 1,
#' the output is arranged in a 3D array.
#' @param seed an object specifying if and how the random number generator should be 
#' initialized (‘seeded’). Either NULL or an integer that will be used in a call to set.seed
#' before simulating the response vectors. If set, the value is saved as the "seed" attribute 
#' of the returned value. The default, NULL will not change the random generator state, and 
#' return .Random.seed as the "seed" attribute, see ‘Value’
#' @param ... additional optional arguments.
#'
#' @details
#' A \code{simulate} function for \code{mlm} objects, which simulates one or more 
#' sets of responses from a Multivariate Linear Model (\code{mlm}) object.
#' If multiple sets of responses are requested, they are returned in a 3D array,
#' with simulation number along the third dimension.
#' 
#' The weights argument is currently ignored -- a constant 
#' variance-covariance matrix assumed for \code{mlm}.
#' 
#' @return A matrix of simulated values for the response (or an array, for \code{nsim} greater than 1)
#' 
#' @author David Warton <david.warton@@unsw.edu.au>
#' 
#' @seealso \code{\link{lm}}, \code{\link{simulate}}
#' @examples
#' # fit a mlm to iris data:
#' data(iris)
#' iris.mlm=lm(cbind(Sepal.Length,Sepal.Width,Petal.Length,Petal.Width)~Species,data=iris)
#' # simulate new responses:
#' simulate(iris.mlm)
#'
#' @method simulate mlm
#' @importFrom stats simulate
#' @export
simulate.mlm = function (object, nsim = 1, seed = NULL, ...) 
{
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
  if (!(is.null(w <- object$weights) || (length(w) == 
                                         1L && w == 1))) vars <- vars/w
  val <- replicate(nsim,ftd+mvtnorm::rmvnorm(n,sigma=vars))
  colnames(val)=nm[[2]]
  rownames(val) <- nm[[1]]
  if(nsim==1)
    val=val[,,1]
  else
    dimnames(val)[[3]] <- paste0("sim_",1:nsim)
  attr(val, "seed") <- RNGstate
  val
}



