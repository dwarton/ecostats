#' Simulate Responses from a (Multivariate) Linear Model
#'
#' Simulate one or more sets of responses from a Linear Model (\code{lm}) or Multivariate Linear Model (\code{mlm}) object.
#'
#' @param object a \code{lm} or \code{mlm} object, typically the result of calling \code{lm} with a matrix response.
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
#' simulate(iris.mlm)
#'
#' @method simulate lm
#' @method simulate mlm
#' @export
simulate.lm = function (object, nsim = 1, seed = NULL, ...) 
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
  fam <- if (isGlm <- inherits(object, "glm")) 
    object$family$family
  else "gaussian"
  ftd <- fitted(object)
  isMlm <- identical(fam, "gaussian") && is.matrix(ftd)
  nm <- if (isMlm) 
    dimnames(ftd)
  else names(ftd)
  if (isMlm) 
  {  # set dimensions of simulated data object
    ftd <- fitted(object)
    n <- nrow(ftd)
    vars <- var(residuals(object))
    #  if (!(is.null(w <- object$weights) || (length(w) == 
    #                                         1L && w == 1))) vars <- vars/w
    val <- replicate(nsim,ftd+mvtnorm::rmvnorm(n,sigma=vars))
    if(nsim==1)
      val=matrix(val,ncol=ncol(ftd))
    else
      dimnames(val)[[3]] = paste0("sim_",seq_len(nsim))
    colnames(val)=nm[[2]]
    rownames(val) = nm[[1]]
  }
  else
  {
    n <- length(ftd)
    ntot <- n * nsim
    val <- switch(fam, gaussian = {
    vars <- deviance(object)/df.residual(object)
    if (isGlm) {
        if (!is.null(object$prior.weights)) vars <- vars/object$prior.weights
      } else if (!(is.null(w <- object$weights) || (length(w) == 
                                                    1L && w == 1))) vars <- vars/w
      ftd + rnorm(ntot, sd = sqrt(vars))
    }
    , if (!is.null(object$family$simulate)) object$family$simulate(object, 
                                                                  nsim) else stop(gettextf("family '%s' not implemented", 
                                                                                           fam), domain = NA))
  if (!is.list(val)) {
    dim(val) <- c(n, nsim)
    val <- as.data.frame(val)
  }
  else class(val) <- "data.frame"
  names(val) <- paste0("sim_", seq_len(nsim))
  if (!is.null(nm)) 
    row.names(val) <- nm
  }
  attr(val, "seed") <- RNGstate
  val
}



#' Simulate Responses
#'
#' Simulate one or more sets of responses from the distribution corresponding to the fitted model object.
#' This function is the same one found in \code{stats} but rewritten here to overcome weird compilation issues.
#'
#' @param object an object representing a fitted model.
#' @param nsim number of replicate datasets to simulate. Defaults to \code{1}.
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
#' @return A vector or matrix of \code{nsim} sets of simulated responses
#' 
#' @author David Warton <david.warton@@unsw.edu.au>
#' 
#' @examples
#' 
#' x <- 1:5
#' mod1 <- lm(c(1:3, 7, 6) ~ x)
#' S1 <- simulate(mod1, nsim = 4)
#' ## repeat the simulation:
#' .Random.seed <- attr(S1, "seed")
#' identical(S1, simulate(mod1, nsim = 4))
#' 
#' S2 <- simulate(mod1, nsim = 200, seed = 101)
#' rowMeans(S2) # should be about the same as
#' fitted(mod1)
#' 
#' ## repeat identically:
#' (sseed <- attr(S2, "seed")) # seed; RNGkind as attribute
#' stopifnot(identical(S2, simulate(mod1, nsim = 200, seed = sseed)))
#' 
#' ## To be sure about the proper RNGkind, e.g., after
#' RNGversion("2.7.0")
#' ## first set the RNG kind, then simulate
#' do.call(RNGkind, attr(sseed, "kind"))
#' identical(S2, simulate(mod1, nsim = 200, seed = sseed))
#' 
#' ## Binomial GLM examples
#' yb1 <- matrix(c(4, 4, 5, 7, 8, 6, 6, 5, 3, 2), ncol = 2)
#' modb1 <- glm(yb1 ~ x, family = binomial)
#' S3 <- simulate(modb1, nsim = 4)
#' # each column of S3 is a two-column matrix.
#' 
#' x2 <- sort(runif(100))
#' yb2 <- rbinom(100, prob = plogis(2*(x2-1)), size = 1)
#' yb2 <- factor(1 + yb2, labels = c("failure", "success"))
#' modb2 <- glm(yb2 ~ x2, family = binomial)
#' S4 <- simulate(modb2, nsim = 4)
#' # each column of S4 is a factor
#' data(iris)
#'
#' @export
simulate = function (object, nsim = 1, seed = NULL, ...) 
  UseMethod("simulate")
