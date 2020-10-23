#' Normal Quantile-Quantile Plots with Simulation Envelopes
#'
#' Produces a normal QQ plot of a \emph{fitted model} \code{y}, with a user-specified 
#' line to compare to "theoretical" quantiles, and global envelopes constructed
#' by simulating new residuals. Global envelopes are constructed using the
#' \code{GET} package for simultaneous control of error rates over the whole curve.
#'
#' @param y can be a set of values for which we wish to check (multivariate) normality, 
#' or it can be \emph{any} object that responds to the \code{residuals} and \code{simulate}
#' functions. The function was designed with models in mind whose residuals would be are
#' approximately normally distributed if the model were correct.
#' @param main the plot title (if a plot is produced)
#' @param xlab \code{x} axis label (if a plot is produced)
#' @param ylab \code{y} axis label (if a plot is produced)
#' @param plot.it logical. Should the result be plotted?
#' @param plot.it logical. Should the result be plotted?
#' @param n.sim the number of simulated sets of residuals to be generated, to which
#'  the observed residuals will be compared. The default is 999 datasets.
#' @param col color of points
#' @param col.envelope the colour of the simulation envelope, defaulting to 
#'  "darkolivegreen". Because it's cool.
#' @param conf.level the confidence level to use in constructing the envelope.
#' @param alpha.f the amount of transparency to use for the shaded region representing
#'  the confidence envelope, values closer to zero are more transparent.
#' @param type the type of global envelope to construct, see 
#'  \code{\link[GET]{global_envelope_test}} for details. Defaults \code{"st"} uses 
#'  studentized envelope tests to protect for unequal variance, which has performed well
#'  in \code{qqnorm} simulations.
#' @param add.line the line to be added to the plot. By default, a one-to-one line is 
#'  added, assuming that residuals are constructed in such a way that they would be
#'  standard normal if assumptions were satisfied. In other cases a better option would be
#'  \code{qqline(y)}. 
#' @param transform a character vector pointing to a function that should be applied to both
#'  theoretical and sample residuals prior constructing an envelope. The main potential
#'  use is to set \code{transform="pnorm"} to map across to PP-plots prior to envelope 
#'  construction. 
#' @param do.fit logical. Should residuals be simulated by constructing new responses
#'  and refitting the model (\code{TRUE}) or should should they be simulated as independent
#'  standard normal random variables (\code{FALSE})? The latter is much faster but is
#'  approximate, and is not appropriate in cases where residuals are not standard normal
#'  when assumptions are satisfied.
#' @param ... further arguments sent through to \code{plot}
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
#' 
#' The best way to use this function, and the default, is to take a model that has
#' been fitted to data and simulate new responses from the model, re-fit the model, and 
#' construct new residuals. This "parametric bootstrap" approach will be computationally 
#' intensive for large or hard-to-fit models. The alternative, \code{do.fit=FALSE}, will
#' simulate independent standard normal data directly and use these to construct the envelope.
#' This faster alternative should only be used if you would in fact expect residuals to be
#' independent and standard normal if model assumptions were satisfied.
#' 
#' The simulated data and subsequent analysis are also used to obtain a \emph{P}-value 
#' for the test that model assumptions are correct. This test uses a 'parametric bootstrap'
#' test (for \code{do.fit=TRUE}), and tests if sample residuals are unusually far from
#' the values expected of them if model assumptions were satisfied. For details see
#' \code{\link[GET]{global_envelope_test}}.#' 
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
#' 
#' @seealso \code{\link{qqnorm}}, \code{\link{qqline}}
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
#' qqenvelope(iris.mlm,n.sim=499)
#' 
#' @export
qqenvelope = function (y, main = "Normal Q-Q Plot", xlab = "Theoretical Quantiles", 
                       ylab = "Sample Quantiles", plot.it = TRUE, 
                       n.sim=999, col=par("col"), col.envelope="darkolivegreen", conf.level=0.95, 
                       alpha.f=0.1, type="st", add.line=c(0,1), transform = NULL,
                       do.fit=TRUE, ...) 
{
  # is y data or model fit? Make y data, object model fit
  if(is.numeric(y))
    object = lm(y~1)
  else
  {
    object=y
    y=model.response(model.frame(object))
  }
  
  # get qq points
  if(inherits(object,"lm")) #use rstudent if an lm object. 
    y=rstudent(object)
  else
    y=residuals(object)
  
  # check if it is multivariate data, if so set a flag for later and vectorise res
  if(length(dim(y))>1)
  {
    nResp=dim(y)[2]
    is.mva = TRUE
  }
  else
  {
    is.mva=FALSE
    nResp=1
  }
  
  qq.x=qqnorm(as.vector(y),plot.it=FALSE)
  n.obs=length(y)
  
  # simulate to get limits around these
  if(do.fit)
  {
    modelF=model.frame(object)
    yNew=simulate(object,n.sim)
    if(is.mva) # for multivariate data, vectorise res for each sim dataset
      yNew=apply(yNew,3,cbind)
    rnorms = matrix(NA,nrow(yNew),ncol(yNew))
    for(iCol in 1:ncol(yNew))
    {
      modelF[1]=matrix(yNew[,iCol],ncol=nResp,dimnames=dimnames(y))
      newFit=try(update(object,data=modelF))
      if(inherits(newFit,"lm")) #use rstudent if an lm object. 
        rnorms[,iCol]=rstudent(newFit)
      else
        rnorms[,iCol]=residuals(newFit)
    }
  }
  else  
    rnorms = matrix(rnorm(n.obs*n.sim),n.obs,n.sim)
  
  qqSort=apply(rnorms,2,sort)
  
  
  # if required, transform data
  if (is.null(transform)==FALSE)
  {
    qq.x$x = do.call(transform,list(qq.x$x))
    qq.x$y = do.call(transform,list(qq.x$y))
    qqSort = do.call(transform,list(qqSort))
    y      = do.call(transform,list(y))
  }
  
  #use the Global Envelope Test package to get global envelope
  ySort=sort(y, index.return=TRUE)
  datCurves=GET::create_curve_set(list(obs=ySort$x,sim_m=qqSort))
  cr=GET::global_envelope_test(datCurves,type=type,alpha=1-conf.level)
  
  if(plot.it==TRUE)
  {
    plot(qq.x$x,qq.x$y,main=main,xlab=xlab,ylab=ylab, col=col, ...)
    lines(range(qq.x$x),add.line[1]+add.line[2]*range(qq.x$x),col=col.envelope)
  }
  
  
  #do a qq plot and add sim envelope
  if(plot.it)
  {
    col_transparent = adjustcolor(col.envelope, alpha.f = alpha.f)
    polygon(sort(qq.x$x)[c(1:n.obs,n.obs:1)],c(cr$lo,cr$hi[n.obs:1]), 
            col=col_transparent, border=NA)
  }
  
  # return a list with qq values and limits, ordered same way as input data
  qqLo = qqHi = rep(NA,n.obs)
  qqLo[ySort$ix] = cr$lo
  qqHi[ySort$ix] = cr$hi
  invisible(list(x = qq.x$x, y = qq.x$y, lo=qqLo, hi=qqHi, p.value=attr(cr,"p")))
}
