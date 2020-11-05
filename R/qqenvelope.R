#' Normal Quantile-Quantile Plots with Global Simulation Envelopes
#'
#' Produces a normal QQ plot of a \emph{fitted model} \code{y}, with a user-specified 
#' line to compare to "theoretical" quantiles, and global envelopes constructed
#' by simulating new residuals. Global envelopes are constructed using the
#' \code{GET} package for simultaneous control of error rates over the whole curve.
#'
#' @param y can be a set of values for which we wish to check (multivariate) normality, 
#' or it can be \emph{any} object that responds to the \code{residuals} and \code{simulate}
#' functions. The function was designed with models in mind whose residuals would be 
#' approximately normally distributed if the model were correct.
#' @param ylab \code{y} axis label (if a plot is produced).
#' @param n.sim the number of simulated sets of residuals to be generated, to which
#'  the observed residuals will be compared. The default is 199 datasets.
#' @param conf.level the confidence level to use in constructing the envelope.
#' @param ... further arguments sent through to \code{plotenvelope}
#' 
#' @details
#' A challenge when interpreting a \code{qqplot} is understanding the extent to which
#' deviations from expected values could be due to random noise (sampling variation)
#' rather than actual assumption violations. This function is intended to assess this, 
#' by simulating multiple realizations of residuals in situations where assumptions 
#' are satisfied, and plotting a simulation envelope around these at level \code{conf.level}.
#' 
#' This function can take data (univariate or multivariate) and check for (multivariate)
#' normality, or it can take a fitted model and use qq plots to interrogate residuals and
#' see if they are behaving as we would expect them to if the model were true.
#' 
#' The envelope is global, constructed using the \code{\link[GET]{GET-package}}, meaning that
#' (for example) a 95% global envelope should contain \emph{all} residuals for 95% of datasets
#' for which assumptions are actually satisfied. So if \emph{any} data points lie outside the
#' envelope we have evidence that assumptions of the fitted model are not satisfied. 
#' The \code{\link[GET]{GET-package}} was originally constructed to place envelopes around functions, motivated by
#' the problem of diagnostic testing of spatial processes (Myllymäki et al 2017), but it can equally
#' well be applied here, by treating the set of residuals (order according to the x-axis) as point-wise evaluations of a function.
#' 
#' For further details refer to \code{\link{plotenvelope}}, which is called to construct the plot.
#' 
#' @return a qqplot with simulation envelope is returned, and additionally:
#' \item{x}{a vector of theoretical quantiles from the standard normal sorted from
#'  smallest to largest}
#' \item{y}{a vector of observed residuals sorted from smallest to largest}
#' \item{lo}{lower bounds on the global simulation envelope for residuals}
#' \item{hi}{upper bounds on the global simulation envelope for residuals}
#' \item{p.value}{A \emph{P}-value for the test that model assumptions are correct,
#'  using a 'parametric bootstrap' test, based on how far sample residuals depart from
#'  the values expected of them if model assumptions were satisfied.}
#' 
#' @author David Warton <david.warton@@unsw.edu.au>
#' @references 
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017), Global envelope tests for spatial processes. J. R. Stat. Soc. B, 79: 381-404. doi:10.1111/rssb.12172
#' 
#' @seealso \code{\link{qqnorm}}, \code{\link{qqline}}, \code{\link{plotenvelope}}
#' @examples
#' # simulate some data and fit a qq plot:
#' y=rnorm(20)
#' qqenvelope(y)
#' 
#' # fit a multivariate linear model to the iris dataset:
#' data(iris)
#' Y = with(iris, cbind(Sepal.Length,Sepal.Width,Petal.Length,Petal.Width))
#' iris.mlm=lm(Y~Species,data=iris)
#' # check normality assumption:
#' qqenvelope(iris.mlm,n.sim=99)
#' 
#' @export
qqenvelope = function (y, n.sim=199, ylab="Sample Quantiles", conf.level=0.95, ...) 
{
  # is y data or model fit? Make y data, object model fit
  if(is.numeric(y))
    y = lm(y~1)
  
  out = plotenvelope(y, which=2, ylab=ylab, n.sim=n.sim, conf.level=conf.level, ...)
  invisible(out[[2]])
}
